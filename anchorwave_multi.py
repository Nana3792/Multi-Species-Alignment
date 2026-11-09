#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path


class PipelineError(RuntimeError):
    pass


@dataclass(frozen=True)
class FastaIndex:
    length: int
    offset: int
    line_bases: int
    line_bytes: int


@dataclass(frozen=True)
class Anchor:
    ref_chr: str
    ref_start: int
    ref_end: int
    query_chr: str
    query_start: int
    query_end: int
    strand: str
    gene: str
    block: str
    score: float | None
    line: int

    @property
    def key(self):
        return self.ref_chr, self.ref_start, self.ref_end


@dataclass(frozen=True)
class SharedAnchor:
    anchor_id: str
    ref_chr: str
    ref_start: int
    ref_end: int
    gene: str
    queries: dict[str, Anchor]


COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhd bn".replace(" ", ""),
)


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def natural_key(text: str):
    return tuple(int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text))


def file_ok(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def full_path(base: Path, value: str) -> Path:
    path = Path(os.path.expandvars(os.path.expanduser(value)))
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def safe_name(name: str, label: str) -> str:
    if not name or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
        raise PipelineError(f"Invalid {label} name: {name!r}")
    return name


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_config(path: Path) -> dict:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise PipelineError(f"Configuration file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise PipelineError(f"Invalid JSON in {path}: line {exc.lineno}, column {exc.colno}") from exc
    if not isinstance(raw, dict):
        raise PipelineError("Configuration root must be a JSON object.")

    path = path.resolve()
    project = full_path(path.parent, str(raw.get("project_dir", ".")))
    ref_raw, query_raw = raw.get("reference"), raw.get("queries")
    params, tools, matching = raw.get("parameters", {}), raw.get("tools", {}), raw.get("anchor_matching", {})
    if not isinstance(ref_raw, dict):
        raise PipelineError("Missing or invalid 'reference' object in config.")
    if not isinstance(query_raw, list) or len(query_raw) < 2:
        raise PipelineError("'queries' must contain at least two query species.")

    ref = {
        "name": safe_name(str(ref_raw.get("name", "reference")), "reference"),
        "fasta": full_path(project, str(ref_raw.get("fasta", ""))),
        "gff3": full_path(project, str(ref_raw.get("gff3", ""))),
    }
    for key in ("cds", "sam"):
        if ref_raw.get(key):
            ref[key] = full_path(project, str(ref_raw[key]))

    queries, names = [], {ref["name"]}
    for i, item in enumerate(query_raw, 1):
        if not isinstance(item, dict):
            raise PipelineError(f"Query entry {i} must be a JSON object.")
        name = safe_name(str(item.get("name", "")), f"query {i}")
        if name in names:
            raise PipelineError(f"Duplicate species name: {name}")
        names.add(name)
        query = {
            "name": name,
            "fasta": full_path(project, str(item.get("fasta", ""))),
            "R": int(item.get("R", params.get("reference_depth", 1))),
            "Q": int(item.get("Q", params.get("query_depth", 1))),
        }
        for key in ("gff3", "sam", "anchors", "maf", "fragments_maf"):
            if item.get(key):
                query[key] = full_path(project, str(item[key]))
        queries.append(query)

    method = str(matching.get("method", "exact")).lower()
    duplicate_policy = str(matching.get("duplicate_policy", "best_score"))
    if method != "exact":
        raise PipelineError("Only anchor_matching.method='exact' is supported.")
    if duplicate_policy != "best_score":
        raise PipelineError("Only duplicate_policy='best_score' is supported.")

    min_exon = params.get("min_exon_length")
    return {
        "config_path": path,
        "project_dir": project,
        "output_dir": full_path(project, str(raw.get("output_dir", "results/awmulti"))),
        "logs_dir": full_path(project, str(raw.get("logs_dir", "logs/awmulti"))),
        "tools": {name: str(tools.get(name, name)) for name in ("anchorwave", "minimap2", "mafft")},
        "reference": ref,
        "queries": queries,
        "parameters": {
            "minimap2_threads": max(1, int(params.get("minimap2_threads", 4))),
            "anchorwave_threads": max(1, int(params.get("anchorwave_threads", 1))),
            "mafft_jobs": max(1, int(params.get("mafft_jobs", 4))),
            "mafft_threads": max(1, int(params.get("mafft_threads", 1))),
            "min_exon_length": None if min_exon is None else max(0, int(min_exon)),
            "include_novel_anchors": bool(params.get("include_novel_anchors", True)),
            "write_block_fastas": bool(params.get("write_block_fastas", True)),
        },
        "anchor_matching": {"method": method, "duplicate_policy": duplicate_policy},
    }


def resolve_tool(spec: str, project: Path) -> str:
    candidate = Path(os.path.expandvars(os.path.expanduser(spec)))
    if candidate.is_absolute() or os.sep in spec or (os.altsep and os.altsep in spec):
        path = candidate if candidate.is_absolute() else project / candidate
        path = path.resolve()
        if not path.is_file() or not os.access(path, os.X_OK):
            raise PipelineError(f"Executable not found or not executable: {path}")
        return str(path)
    path = shutil.which(spec)
    if not path:
        raise PipelineError(f"Executable not found on PATH: {spec}")
    return path


def tool_version(exe: str, *options: list[str]) -> str:
    for option in options:
        try:
            result = subprocess.run([exe, *option], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, timeout=20, check=False)
        except (OSError, subprocess.TimeoutExpired):
            continue
        if result.stdout.strip():
            return result.stdout.strip().splitlines()[0][:500]
    return "unknown"


def fasta_names(path: Path) -> set[str]:
    names = set()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                if not name or name in names:
                    raise PipelineError(f"Invalid or duplicate FASTA header in {path}: {name!r}")
                names.add(name)
    if not names:
        raise PipelineError(f"No FASTA records found in {path}")
    return names


def gff_names(path: Path) -> set[str]:
    names = set()
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                raise PipelineError(f"Malformed GFF3 line {line_no} in {path}")
            names.add(fields[0])
    if not names:
        raise PipelineError(f"No feature sequence IDs found in GFF3: {path}")
    return names


def validate_inputs(config: dict) -> dict:
    files = [("reference FASTA", config["reference"]["fasta"]),
             ("reference GFF3", config["reference"]["gff3"])]
    for query in config["queries"]:
        files.append((f"{query['name']} FASTA", query["fasta"]))
        if query.get("gff3"):
            files.append((f"{query['name']} GFF3", query["gff3"]))
    for label, path in files:
        if not file_ok(path):
            raise PipelineError(f"Missing or empty {label}: {path}")
        if path.suffix == ".gz":
            raise PipelineError(f"Compressed input is not supported; decompress {label}: {path}")

    config["resolved_tools"] = {name: resolve_tool(spec, config["project_dir"])
                                for name, spec in config["tools"].items()}
    missing = sorted(gff_names(config["reference"]["gff3"]) - fasta_names(config["reference"]["fasta"]),
                     key=natural_key)
    if missing:
        preview = ", ".join(missing[:10])
        suffix = "" if len(missing) <= 10 else f" ... and {len(missing) - 10} more"
        raise PipelineError(f"Reference: {len(missing)} GFF3 sequence IDs are missing from the FASTA: {preview}{suffix}")

    for query in config["queries"]:
        fasta = fasta_names(query["fasta"])
        if query.get("gff3"):
            missing = gff_names(query["gff3"]) - fasta
            if missing:
                raise PipelineError(f"{query['name']}: {len(missing)} GFF3 sequence IDs are missing from the FASTA")

    tools = config["resolved_tools"]
    return {
        "anchorwave": tool_version(tools["anchorwave"], ["--version"], ["-h"]),
        "minimap2": tool_version(tools["minimap2"], ["--version"]),
        "mafft": tool_version(tools["mafft"], ["--version"]),
        "python": sys.version.splitlines()[0],
    }


def run_command(command: list[str], log: Path, output: Path | None = None, dry_run: bool = False) -> float:
    log.parent.mkdir(parents=True, exist_ok=True)
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
    text = " ".join(shlex.quote(str(x)) for x in command)
    if dry_run:
        print("DRY RUN:", text)
        return 0.0

    start = time.perf_counter()
    with log.open("wt", encoding="utf-8") as err:
        err.write(f"# command: {text}\n")
        if output:
            with output.open("wb") as out:
                result = subprocess.run(command, stdout=out, stderr=err, check=False)
        else:
            result = subprocess.run(command, stdout=err, stderr=subprocess.STDOUT, check=False)
        elapsed = time.perf_counter() - start
        err.write(f"\n# elapsed_seconds: {elapsed:.3f}\n# return_code: {result.returncode}\n")
    if result.returncode:
        raise PipelineError(f"Command failed with exit code {result.returncode}: {text}\nLog: {log}")
    return elapsed


class IndexedFasta:
    def __init__(self, path: Path):
        self.path = path
        self.index = self._index()

    def _index(self) -> dict[str, FastaIndex]:
        cache, stat = self.path.with_name(self.path.name + ".awmulti-index.json"), self.path.stat()
        if cache.is_file():
            try:
                data = json.loads(cache.read_text(encoding="utf-8"))
                if data.get("size") == stat.st_size and data.get("mtime_ns") == stat.st_mtime_ns:
                    return {name: FastaIndex(**values) for name, values in data["records"].items()}
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                pass

        records = {}
        with self.path.open("rb") as handle:
            name = None
            length = offset = line_bases = line_bytes = 0

            def finish():
                if name is not None:
                    if not length or not line_bases or not line_bytes:
                        raise PipelineError(f"FASTA record {name!r} has no sequence in {self.path}")
                    records[name] = FastaIndex(length, offset, line_bases, line_bytes)

            while True:
                line_start, line = handle.tell(), handle.readline()
                if not line:
                    finish(); break
                if line.startswith(b">"):
                    finish()
                    header = line[1:].strip().split(None, 1)[0]
                    if not header:
                        raise PipelineError(f"Empty FASTA header in {self.path}")
                    name = header.decode("utf-8")
                    if name in records:
                        raise PipelineError(f"Duplicate FASTA record {name!r} in {self.path}")
                    length = offset = line_bases = line_bytes = 0
                elif name is None:
                    if line.strip():
                        raise PipelineError(f"Sequence data before first FASTA header in {self.path}")
                else:
                    seq = line.rstrip(b"\r\n")
                    if seq:
                        if not line_bases:
                            line_bases, line_bytes, offset = len(seq), len(line), line_start
                        length += len(seq)
        if not records:
            raise PipelineError(f"No FASTA records found in {self.path}")

        cache.write_text(json.dumps({
            "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
            "records": {name: r.__dict__ for name, r in records.items()},
        }, indent=2, sort_keys=True), encoding="utf-8")
        return records

    def fetch(self, name: str, start: int, end: int) -> str:
        record = self.index.get(name)
        if record is None:
            raise PipelineError(f"Sequence {name!r} not found in FASTA {self.path}")
        if start < 1 or end < start or end > record.length:
            raise PipelineError(f"Invalid interval {name}:{start}-{end} in {self.path}")

        current, remaining, seq = start - 1, end - start + 1, bytearray()
        with self.path.open("rb") as handle:
            while remaining:
                line, col = divmod(current, record.line_bases)
                take = min(remaining, record.line_bases - col)
                handle.seek(record.offset + line * record.line_bytes + col)
                chunk = handle.read(take)
                if len(chunk) != take:
                    raise PipelineError(f"Unexpected end of FASTA while reading {name}:{start}-{end}")
                seq.extend(chunk); current += take; remaining -= take
        return seq.decode("ascii").upper()


def parse_anchor_file(path: Path, include_novel: bool = True) -> list[Anchor]:
    anchors, header = [], False
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if fields[0] == "refChr":
                header = True; continue
            if len(fields) < 10:
                raise PipelineError(f"Malformed AnchorWave row in {path} at line {line_no}")
            if fields[7] == "interanchor" or (not include_novel and "localAlignment" in fields[7]):
                continue
            try:
                anchor = Anchor(fields[0], int(fields[1]), int(fields[2]), fields[3], int(fields[4]), int(fields[5]),
                                fields[6], fields[7], fields[8].strip() or "NA",
                                None if fields[9] == "NA" else float(fields[9]), line_no)
            except ValueError as exc:
                raise PipelineError(f"Invalid numeric value in {path} at line {line_no}: {text}") from exc
            if anchor.strand not in {"+", "-"}:
                raise PipelineError(f"Invalid strand {anchor.strand!r} in {path} at line {line_no}")
            if anchor.ref_start < 1 or anchor.ref_end < anchor.ref_start or anchor.query_start < 1 or anchor.query_end < anchor.query_start:
                raise PipelineError(f"Invalid anchor interval in {path} at line {line_no}")
            anchors.append(anchor)
    if not header or not anchors:
        raise PipelineError(f"No usable AnchorWave anchors found in {path}")
    return anchors


def anchor_rank(a: Anchor):
    return (float("-inf") if a.score is None else a.score,
            a.query_end - a.query_start + 1, a.query_chr, -a.query_start, -a.query_end, -a.line)


def select_anchors(records: list[Anchor]):
    groups = {}
    for record in records:
        groups.setdefault(record.key, []).append(record)
    return ({key: max(values, key=anchor_rank) for key, values in groups.items()},
            {key: values for key, values in groups.items() if len(values) > 1})


def build_shared_anchors(query_records: dict[str, list[Anchor]]):
    selected, duplicates = {}, {}
    for name, records in query_records.items():
        selected[name], duplicates[name] = select_anchors(records)
    keys = set.intersection(*(set(x) for x in selected.values()))
    shared = []
    for i, key in enumerate(sorted(keys, key=lambda x: (natural_key(x[0]), x[1], x[2])), 1):
        records = {name: selected[name][key] for name in selected}
        genes = [record.gene for record in records.values()]
        gene = max(set(genes), key=lambda value: (genes.count(value), value))
        shared.append(SharedAnchor(f"anchor_{i:06d}", *key, gene, records))
    return shared, duplicates


def write_shared_table(path: Path, shared: list[SharedAnchor], queries: list[str]) -> None:
    header = ["anchor_id", "ref_chr", "ref_start", "ref_end", "gene"]
    for name in queries:
        header += [f"{name}_chr", f"{name}_start", f"{name}_end", f"{name}_strand", f"{name}_score", f"{name}_block"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(header)
        for anchor in shared:
            row = [anchor.anchor_id, anchor.ref_chr, anchor.ref_start, anchor.ref_end, anchor.gene]
            for name in queries:
                r = anchor.queries[name]
                row += [r.query_chr, r.query_start, r.query_end, r.strand,
                        "NA" if r.score is None else f"{r.score:.10g}", r.block]
            writer.writerow(row)


def write_duplicate_table(path: Path, duplicates: dict) -> None:
    header = ["query", "ref_chr", "ref_start", "ref_end", "candidate_count", "selected_query_chr",
              "selected_query_start", "selected_query_end", "selected_strand", "selected_score", "all_candidates"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(header)
        for query in sorted(duplicates):
            for key in sorted(duplicates[query], key=lambda x: (natural_key(x[0]), x[1], x[2])):
                values, best = duplicates[query][key], max(duplicates[query][key], key=anchor_rank)
                all_candidates = ";".join(f"{r.query_chr}:{r.query_start}-{r.query_end}:{r.strand}:{'NA' if r.score is None else r.score}" for r in values)
                writer.writerow([query, *key, len(values), best.query_chr, best.query_start, best.query_end,
                                 best.strand, "NA" if best.score is None else best.score, all_candidates])


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="ascii") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            handle.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))


def read_fasta(path: Path) -> dict[str, str]:
    records, name = {}, None
    with path.open("rt", encoding="ascii") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                name = text[1:].split()[0]
                if name in records:
                    raise PipelineError(f"Duplicate FASTA name {name!r} in {path}")
                records[name] = []
            elif name is None:
                raise PipelineError(f"Sequence before header in {path} at line {line_no}")
            else:
                records[name].append(text)
    if not records:
        raise PipelineError(f"No records in FASTA {path}")
    return {name: "".join(parts) for name, parts in records.items()}


def extract_blocks(shared: list[SharedAnchor], reference_name: str, reference_fasta: Path,
                   query_fastas: dict[str, Path], directory: Path, force: bool) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    ref = IndexedFasta(reference_fasta)
    queries = {name: IndexedFasta(path) for name, path in query_fastas.items()}
    for anchor in shared:
        out = directory / f"{anchor.anchor_id}.fa"
        if file_ok(out) and not force:
            continue
        records = [(reference_name, ref.fetch(anchor.ref_chr, anchor.ref_start, anchor.ref_end))]
        for name, fasta in queries.items():
            r = anchor.queries[name]
            seq = fasta.fetch(r.query_chr, r.query_start, r.query_end)
            records.append((name, reverse_complement(seq) if r.strand == "-" else seq))
        write_fasta(out, records)


def run_mafft(mafft: str, shared: list[SharedAnchor], blocks: Path, aligned: Path, logs: Path,
              jobs: int, threads: int, force: bool) -> float:
    aligned.mkdir(parents=True, exist_ok=True); logs.mkdir(parents=True, exist_ok=True)

    def one(anchor: SharedAnchor) -> float:
        inp = blocks / f"{anchor.anchor_id}.fa"
        out = aligned / f"{anchor.anchor_id}.aligned.fa"
        if file_ok(out) and not force:
            return 0.0
        return run_command([mafft, "--auto", "--thread", str(threads), str(inp)],
                           logs / f"{anchor.anchor_id}.mafft.log", out)

    total = 0.0
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(one, anchor): anchor.anchor_id for anchor in shared}
        for i, future in enumerate(as_completed(futures), 1):
            try:
                total += future.result()
            except Exception as exc:
                raise PipelineError(f"MAFFT failed for {futures[future]}: {exc}") from exc
            if i % 100 == 0 or i == len(shared):
                print(f"MAFFT blocks completed: {i}/{len(shared)}")
    return total


def concatenate(shared: list[SharedAnchor], species: list[str], aligned: Path,
                output: Path, partitions: Path, force: bool) -> None:
    if file_ok(output) and file_ok(partitions) and not force:
        return
    output.parent.mkdir(parents=True, exist_ok=True); partitions.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="awmulti_concat_") as tmp:
        temp = {name: Path(tmp) / f"{name}.raw" for name in species}
        handles = {name: path.open("wt", encoding="ascii") for name, path in temp.items()}
        try:
            position = 1
            with partitions.open("wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["anchor_id", "ref_chr", "ref_start", "ref_end", "alignment_start", "alignment_end", "alignment_length"])
                for anchor in shared:
                    path = aligned / f"{anchor.anchor_id}.aligned.fa"
                    records = read_fasta(path)
                    if set(records) != set(species):
                        raise PipelineError(f"Unexpected species in {path}")
                    lengths = {len(records[name]) for name in species}
                    if len(lengths) != 1:
                        raise PipelineError(f"Aligned sequences have different lengths in {path}")
                    length = next(iter(lengths))
                    for name in species:
                        handles[name].write(records[name])
                    writer.writerow([anchor.anchor_id, anchor.ref_chr, anchor.ref_start, anchor.ref_end,
                                     position, position + length - 1, length])
                    position += length
        finally:
            for handle in handles.values():
                handle.close()

        with output.open("wt", encoding="ascii") as out:
            for name in species:
                out.write(f">{name}\n")
                with temp[name].open("rt", encoding="ascii") as src:
                    buffer = ""
                    for chunk in iter(lambda: src.read(1024 * 1024), ""):
                        buffer += chunk
                        while len(buffer) >= 80:
                            out.write(buffer[:80] + "\n"); buffer = buffer[80:]
                    if buffer:
                        out.write(buffer + "\n")


def stage_paths(config: dict) -> dict[str, Path]:
    out, logs = config["output_dir"], config["logs_dir"]
    return {
        "reference": out / "reference", "pairwise": out / "pairwise", "multi": out / "multispecies",
        "log_reference": logs / "reference", "log_pairwise": logs / "pairwise", "log_mafft": logs / "mafft",
    }


def reference_outputs(config: dict) -> dict[str, Path]:
    directory = stage_paths(config)["reference"]
    ref = config["reference"]
    return {"cds": ref.get("cds", directory / "reference_cds.fa"), "sam": ref.get("sam", directory / "reference.sam")}


def query_outputs(config: dict, name: str) -> dict[str, Path]:
    p = stage_paths(config); query = next(q for q in config["queries"] if q["name"] == name)
    out, log = p["pairwise"] / name, p["log_pairwise"] / name
    return {
        "dir": out, "log_dir": log, "sam": query.get("sam", out / "query.sam"),
        "anchors": query.get("anchors", out / "anchors.tsv"), "maf": query.get("maf", out / "alignment.maf"),
        "fragments_maf": query.get("fragments_maf", out / "fragments.maf"),
        "minimap_log": log / "minimap2.log", "proali_log": log / "proali.log",
    }


def minimap_command(exe: str, genome: Path, cds: Path, threads: int) -> list[str]:
    return [exe, "-x", "splice", "-t", str(threads), "-k", "12", "-a", "-p", "0.4", "-N", "20", str(genome), str(cds)]


def run_pairwise(config: dict, force: bool, dry_run: bool) -> list[dict]:
    tools, params, p = config["resolved_tools"], config["parameters"], stage_paths(config)
    p["reference"].mkdir(parents=True, exist_ok=True); p["log_reference"].mkdir(parents=True, exist_ok=True)
    ref = reference_outputs(config); cds, ref_sam = ref["cds"], ref["sam"]
    cds.parent.mkdir(parents=True, exist_ok=True); ref_sam.parent.mkdir(parents=True, exist_ok=True)
    timings = []

    def record(stage: str, elapsed: float, status: str):
        timings.append({"stage": stage, "seconds": elapsed, "status": status})

    if force or not file_ok(cds):
        cmd = [tools["anchorwave"], "gff2seq", "-i", str(config["reference"]["gff3"]),
               "-r", str(config["reference"]["fasta"]), "-o", str(cds)]
        if params["min_exon_length"] is not None:
            cmd += ["-m", str(params["min_exon_length"])]
        record("gff2seq", run_command(cmd, p["log_reference"] / "gff2seq.log", dry_run=dry_run), "ran")
    else:
        record("gff2seq", 0.0, "skipped")

    if force or not file_ok(ref_sam):
        record("map_reference", run_command(minimap_command(tools["minimap2"], config["reference"]["fasta"], cds,
                                                             params["minimap2_threads"]),
                                            p["log_reference"] / "minimap2_reference.log", ref_sam, dry_run), "ran")
    else:
        record("map_reference", 0.0, "skipped")

    for query in config["queries"]:
        out = query_outputs(config, query["name"])
        out["dir"].mkdir(parents=True, exist_ok=True); out["log_dir"].mkdir(parents=True, exist_ok=True)
        for key in ("sam", "anchors", "maf", "fragments_maf"):
            out[key].parent.mkdir(parents=True, exist_ok=True)

        stage = f"map_{query['name']}"
        if force or not file_ok(out["sam"]):
            record(stage, run_command(minimap_command(tools["minimap2"], query["fasta"], cds, params["minimap2_threads"]),
                                      out["minimap_log"], out["sam"], dry_run), "ran")
        else:
            record(stage, 0.0, "skipped")

        stage = f"proali_{query['name']}"
        complete = all(file_ok(out[key]) for key in ("anchors", "maf", "fragments_maf"))
        if force or not complete:
            if not dry_run:
                for key in ("anchors", "maf", "fragments_maf"):
                    out[key].unlink(missing_ok=True)
            cmd = [
                tools["anchorwave"], "proali", "-i", str(config["reference"]["gff3"]), "-as", str(cds),
                "-r", str(config["reference"]["fasta"]), "-a", str(out["sam"]), "-ar", str(ref_sam),
                "-s", str(query["fasta"]), "-n", str(out["anchors"]), "-R", str(query["R"]), "-Q", str(query["Q"]),
                "-o", str(out["maf"]), "-f", str(out["fragments_maf"]), "-t", str(params["anchorwave_threads"]),
            ]
            if params["min_exon_length"] is not None:
                cmd += ["-m", str(params["min_exon_length"])]
            record(stage, run_command(cmd, out["proali_log"], dry_run=dry_run), "ran")
        else:
            record(stage, 0.0, "skipped")
    return timings


def run_merge(config: dict, force: bool, dry_run: bool) -> tuple[list[dict], dict]:
    if dry_run:
        print("DRY RUN: parse anchors, intersect coordinates, extract blocks, run MAFFT, concatenate")
        return [{"stage": "merge_and_msa", "seconds": 0.0, "status": "dry-run"}], {}

    params, p = config["parameters"], stage_paths(config); multi = p["multi"]
    multi.mkdir(parents=True, exist_ok=True)
    query_names = [q["name"] for q in config["queries"]]
    records, counts = {}, {}
    for query in config["queries"]:
        path = query_outputs(config, query["name"])["anchors"]
        if not file_ok(path):
            raise PipelineError(f"Missing pairwise anchors for {query['name']}: {path}")
        records[query["name"]] = parse_anchor_file(path, params["include_novel_anchors"])
        counts[query["name"]] = len(records[query["name"]])

    start = time.perf_counter(); shared, duplicates = build_shared_anchors(records); merge_time = time.perf_counter() - start
    if not shared:
        raise PipelineError("No exact shared reference anchors were found across all query species.")

    shared_table, duplicate_table = multi / "shared_anchors_exact.tsv", multi / "duplicate_anchor_candidates.tsv"
    write_shared_table(shared_table, shared, query_names); write_duplicate_table(duplicate_table, duplicates)

    blocks = multi / "blocks"; start = time.perf_counter()
    extract_blocks(shared, config["reference"]["name"], config["reference"]["fasta"],
                   {q["name"]: q["fasta"] for q in config["queries"]}, blocks, force)
    extract_time = time.perf_counter() - start

    aligned = multi / "aligned_blocks"
    mafft_time = run_mafft(config["resolved_tools"]["mafft"], shared, blocks, aligned, p["log_mafft"],
                           params["mafft_jobs"], params["mafft_threads"], force)

    final, partitions = multi / "all_shared_anchors.aligned.fa", multi / "anchor_partitions.tsv"
    start = time.perf_counter(); concatenate(shared, [config["reference"]["name"], *query_names], aligned, final, partitions, force)
    concat_time = time.perf_counter() - start
    if not params["write_block_fastas"]:
        shutil.rmtree(blocks, ignore_errors=True)

    summary = {
        "pairwise_anchor_counts_excluding_interanchors": counts,
        "shared_anchor_count": len(shared),
        "duplicate_reference_keys_by_query": {name: len(values) for name, values in duplicates.items()},
        "outputs": {"shared_anchor_table": str(shared_table), "duplicate_report": str(duplicate_table),
                    "aligned_blocks": str(aligned), "supermatrix": str(final), "partitions": str(partitions)},
    }
    timings = [
        {"stage": "intersect_shared_anchors", "seconds": merge_time, "status": "ran"},
        {"stage": "extract_anchor_sequences", "seconds": extract_time, "status": "ran"},
        {"stage": "mafft_blocks_total_process_time", "seconds": mafft_time, "status": "ran"},
        {"stage": "concatenate_alignments", "seconds": concat_time, "status": "ran"},
    ]
    return timings, summary


def print_plan(config: dict) -> None:
    ref = reference_outputs(config)
    print("Project:", config["project_dir"]); print("Reference:", config["reference"]["name"])
    print("  FASTA:", config["reference"]["fasta"]); print("  GFF3: ", config["reference"]["gff3"])
    print("  CDS:  ", ref["cds"]); print("  SAM:  ", ref["sam"]); print("Queries:")
    for query in config["queries"]:
        out = query_outputs(config, query["name"])
        print(f"  {query['name']}: {query['fasta']}  (-R {query['R']} -Q {query['Q']})")
        print(f"    SAM: {out['sam']}\n    anchors: {out['anchors']}")
    print("Output:", config["output_dir"]); print("Logs:  ", config["logs_dir"])
    print("Final alignment:", stage_paths(config)["multi"] / "all_shared_anchors.aligned.fa")


def write_manifest(config: dict, versions: dict, timings: list[dict], summary: dict) -> Path:
    ref = reference_outputs(config)
    def qrow(q):
        out = query_outputs(config, q["name"])
        return {"name": q["name"], "fasta": str(q["fasta"]), "gff3": str(q["gff3"]) if q.get("gff3") else None,
                "R": q["R"], "Q": q["Q"], "sam": str(out["sam"]), "anchors": str(out["anchors"]),
                "maf": str(out["maf"]), "fragments_maf": str(out["fragments_maf"])}
    manifest = {
        "created_at_local": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "config_file": str(config["config_path"]),
        "config_sha256": sha256(config["config_path"]), "versions": versions,
        "reference": {"name": config["reference"]["name"], "fasta": str(config["reference"]["fasta"]),
                      "gff3": str(config["reference"]["gff3"]), "cds": str(ref["cds"]), "sam": str(ref["sam"])},
        "queries": [qrow(q) for q in config["queries"]], "parameters": config["parameters"],
        "anchor_matching": config["anchor_matching"], "timings": timings, "summary": summary,
    }
    path = config["output_dir"] / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a reference-centred multi-species alignment from AnchorWave results.")
    parser.add_argument("command", choices=["check", "pairwise", "merge", "run"])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        config = read_config(args.config); print_plan(config); versions = validate_inputs(config)
        print("Tool versions:"); [print(f"  {name}: {value}") for name, value in versions.items()]
        if args.command == "check":
            print("Input and tool checks passed."); return 0

        config["output_dir"].mkdir(parents=True, exist_ok=True); config["logs_dir"].mkdir(parents=True, exist_ok=True)
        timings, summary = [], {}
        if args.command in {"pairwise", "run"}:
            timings += run_pairwise(config, args.force, args.dry_run)
        if args.command in {"merge", "run"}:
            more, summary = run_merge(config, args.force, args.dry_run); timings += more
        if not args.dry_run:
            print("Manifest:", write_manifest(config, versions, timings, summary))
        if summary:
            print("Shared anchors:", summary["shared_anchor_count"])
            print("Multiple alignment:", summary["outputs"]["supermatrix"])
            print("Partitions:", summary["outputs"]["partitions"])
        print("Completed successfully."); return 0
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr); return 2
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr); return 130


if __name__ == "__main__":
    raise SystemExit(main())
