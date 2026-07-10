#!/usr/bin/env python3
"""Fast sequence-only contact screening for REBASEfinder motif pairs.

This is a triage step, not a coordinate predictor. It uses the contact head of
an ESM-2 language model to estimate whether two role-specific motif blocks are
likely to contact. Strong results can prioritize candidates for ESMFold or
another coordinate-producing model, but must not be reported as an Angstrom
distance or direct 3D verification.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
from collections import defaultdict
from pathlib import Path


PAIR_RULES = (
    ("amino_mtase_I-IV", "SAM", "Amino-IV"),
    ("c5_mtase_PC-ENV", "N5C-PC", "N5C-ENV"),
    ("mmei_mtase_I-IV", "MmeI-I", "Amino-IV"),
    ("typeI_motor_WA-WB", "P-loop", "HsdR-WB"),
    ("typeI_motor_WB-MIII", "HsdR-WB", "HsdR-MIII"),
    ("typeIII_motor_WA-WB", "ResIII-WA", "ResIII-WB"),
    ("typeIII_motor_WB-MIII", "ResIII-WB", "ResIII-MIII"),
    ("pld_HKD-pair", "PLD-HKD", "PLD-HKD"),
)

OUTPUT_FIELDS = (
    "locus_tag",
    "family_id",
    "enzyme_role",
    "pair_id",
    "motif_a",
    "motif_b",
    "motif_a_start",
    "motif_a_end",
    "motif_b_start",
    "motif_b_end",
    "sequence_length",
    "crop_start",
    "crop_end",
    "crop_length",
    "model",
    "device",
    "max_contact_probability",
    "mean_contact_probability",
    "top3_mean_contact_probability",
    "contact_pairs_ge_0_5",
    "separation_matched_percentile",
    "sequence_separation",
    "contact_status",
    "contact_reason",
    "inference_seconds",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Estimate contacts between REBASEfinder motif pairs directly from "
            "protein sequence with ESM-2. Results are screening hints only."
        )
    )
    parser.add_argument("--fasta", required=True, help="REBASEfinder candidate protein FASTA")
    parser.add_argument("--motifs", required=True, help="DNMB_REBASEfinder_motif_hits.tsv")
    parser.add_argument("--out", required=True, help="Output TSV")
    parser.add_argument(
        "--model",
        default="esm2_t6_8M_UR50D",
        help=(
            "fair-esm pretrained factory name. The 8M default is the fastest; "
            "esm2_t30_150M_UR50D is the recommended confirmation model."
        ),
    )
    parser.add_argument("--device", choices=("auto", "cuda", "mps", "cpu"), default="auto")
    parser.add_argument("--max-aa", type=int, default=1022, help="Maximum ESM-2 sequence length")
    parser.add_argument("--flank", type=int, default=80, help="Flank used when a long sequence must be cropped")
    parser.add_argument("--limit", type=int, default=0, help="Maximum genes to evaluate; 0 means all")
    return parser.parse_args()


def read_fasta(path: str) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    key = None
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]
                records.setdefault(key, [])
            elif key is not None:
                records[key].append(line)
    return {
        key: re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", "".join(parts).upper())
        for key, parts in records.items()
    }


def as_int(value: str | None) -> int | None:
    try:
        number = int(float(value or ""))
    except (TypeError, ValueError):
        return None
    return number if number >= 1 else None


def read_motifs(path: str) -> dict[str, list[dict[str, str]]]:
    by_locus: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(path, encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus = (row.get("locus_tag") or "").strip()
            start = as_int(row.get("start_aa"))
            end = as_int(row.get("end_aa"))
            evidence = (row.get("evidence_level") or "").strip()
            if not locus or start is None or end is None or end < start:
                continue
            if evidence and evidence not in {"supported", "sequence_hint"}:
                continue
            row["_start"] = str(start)
            row["_end"] = str(end)
            by_locus[locus].append(row)
    return by_locus


def motif_pairs(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    pairs: list[dict[str, object]] = []
    for pair_id, motif_a, motif_b in PAIR_RULES:
        hits_a = [row for row in rows if row.get("motif") == motif_a]
        hits_b = [row for row in rows if row.get("motif") == motif_b]
        if not hits_a or not hits_b:
            continue
        for index_a, hit_a in enumerate(hits_a):
            for index_b, hit_b in enumerate(hits_b):
                if motif_a == motif_b and index_a >= index_b:
                    continue
                pairs.append(
                    {
                        "pair_id": pair_id,
                        "motif_a": motif_a,
                        "motif_b": motif_b,
                        "hit_a": hit_a,
                        "hit_b": hit_b,
                    }
                )
    return pairs


def crop_bounds(sequence_length: int, pairs: list[dict[str, object]], max_aa: int, flank: int) -> tuple[int, int] | None:
    if sequence_length <= max_aa:
        return 1, sequence_length
    starts = []
    ends = []
    for pair in pairs:
        for key in ("hit_a", "hit_b"):
            hit = pair[key]
            starts.append(int(hit["_start"]))
            ends.append(int(hit["_end"]))
    motif_start = min(starts)
    motif_end = max(ends)
    if motif_end - motif_start + 1 > max_aa:
        return None
    start = max(1, motif_start - flank)
    end = min(sequence_length, motif_end + flank)
    if end - start + 1 > max_aa:
        extra = end - start + 1 - max_aa
        trim_left = min(extra // 2, motif_start - start)
        start += trim_left
        end -= extra - trim_left
    if end - start + 1 < max_aa:
        missing = max_aa - (end - start + 1)
        add_left = min(missing // 2, start - 1)
        start -= add_left
        end = min(sequence_length, end + missing - add_left)
        start = max(1, end - max_aa + 1)
    return start, end


def choose_device(torch, requested: str):
    if requested == "cuda":
        if not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is unavailable")
        return torch.device("cuda")
    if requested == "mps":
        if not torch.backends.mps.is_available():
            raise RuntimeError("MPS was requested but is unavailable")
        return torch.device("mps")
    if requested == "cpu":
        return torch.device("cpu")
    if torch.cuda.is_available():
        return torch.device("cuda")
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def contact_metrics(np, contacts, a_start: int, a_end: int, b_start: int, b_end: int) -> dict[str, float | int]:
    block = contacts[a_start - 1 : a_end, b_start - 1 : b_end]
    values = block.reshape(-1)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise RuntimeError("motif contact block had no finite values")
    ordered = np.sort(finite)
    top3 = ordered[-min(3, ordered.size) :]
    center_a = (a_start + a_end) / 2.0
    center_b = (b_start + b_end) / 2.0
    separation = abs(center_b - center_a)
    ii, jj = np.triu_indices(contacts.shape[0], k=24)
    separations = jj - ii
    tolerance = max(12.0, separation * 0.25)
    matched = np.abs(separations - separation) <= tolerance
    background = contacts[ii[matched], jj[matched]]
    background = background[np.isfinite(background)]
    if background.size < 100:
        background = contacts[ii, jj]
        background = background[np.isfinite(background)]
    maximum = float(np.max(finite))
    percentile = float(100.0 * np.mean(background <= maximum)) if background.size else math.nan
    return {
        "max": maximum,
        "mean": float(np.mean(finite)),
        "top3": float(np.mean(top3)),
        "n_ge_0_5": int(np.sum(finite >= 0.5)),
        "percentile": percentile,
        "separation": int(round(separation)),
    }


def classify_contact(metrics: dict[str, float | int]) -> tuple[str, str]:
    maximum = float(metrics["max"])
    top3 = float(metrics["top3"])
    percentile = float(metrics["percentile"])
    if maximum >= 0.50 and top3 >= 0.20 and percentile >= 99.0:
        return "contact_strong", "max>=0.50; top3_mean>=0.20; matched_percentile>=99"
    if maximum >= 0.30 and top3 >= 0.12 and percentile >= 95.0:
        return "contact_supportive", "max>=0.30; top3_mean>=0.12; matched_percentile>=95"
    return "contact_weak", "sequence-only contact thresholds not met"


def blank_result(locus: str, pair: dict[str, object], sequence_length: int, model: str, reason: str) -> dict[str, object]:
    hit_a = pair["hit_a"]
    hit_b = pair["hit_b"]
    return {
        "locus_tag": locus,
        "family_id": hit_a.get("family_id", ""),
        "enzyme_role": hit_a.get("enzyme_role", ""),
        "pair_id": pair["pair_id"],
        "motif_a": pair["motif_a"],
        "motif_b": pair["motif_b"],
        "motif_a_start": hit_a["_start"],
        "motif_a_end": hit_a["_end"],
        "motif_b_start": hit_b["_start"],
        "motif_b_end": hit_b["_end"],
        "sequence_length": sequence_length,
        "model": model,
        "contact_status": "contact_unavailable",
        "contact_reason": reason,
    }


def main() -> int:
    args = parse_args()
    sequences = read_fasta(args.fasta)
    motifs = read_motifs(args.motifs)
    work = []
    for locus, rows in motifs.items():
        pairs = motif_pairs(rows)
        if pairs and locus in sequences:
            work.append((locus, sequences[locus], pairs))
    work.sort(key=lambda item: item[0])
    if args.limit > 0:
        work = work[: args.limit]

    try:
        import numpy as np
        import torch
        import esm
    except ImportError as exc:
        raise SystemExit(
            "ESM-2 dependencies are unavailable. Run this script with a Python "
            "environment containing torch and fair-esm (import name: esm). "
            f"Original error: {exc}"
        )

    factory = getattr(esm.pretrained, args.model, None)
    if factory is None or not callable(factory):
        raise SystemExit(f"Unknown fair-esm pretrained model factory: {args.model}")
    device = choose_device(torch, args.device)
    model, alphabet = factory()
    model = model.eval().to(device)
    batch_converter = alphabet.get_batch_converter()
    results: list[dict[str, object]] = []

    for locus, sequence, pairs in work:
        bounds = crop_bounds(len(sequence), pairs, args.max_aa, args.flank)
        if bounds is None:
            for pair in pairs:
                results.append(
                    blank_result(locus, pair, len(sequence), args.model, "motif span exceeds --max-aa")
                )
            continue
        crop_start, crop_end = bounds
        crop = sequence[crop_start - 1 : crop_end]
        inference_started = time.perf_counter()
        try:
            _, _, tokens = batch_converter([(locus, crop)])
            tokens = tokens.to(device)
            with torch.no_grad():
                output = model(tokens, return_contacts=True)
            if device.type == "cuda":
                torch.cuda.synchronize()
            elif device.type == "mps":
                torch.mps.synchronize()
            contacts = output["contacts"][0, : len(crop), : len(crop)].float().cpu().numpy()
            inference_seconds = time.perf_counter() - inference_started
        except Exception as exc:
            if device.type != "cpu":
                device = torch.device("cpu")
                model = model.to(device)
                _, _, tokens = batch_converter([(locus, crop)])
                with torch.no_grad():
                    output = model(tokens, return_contacts=True)
                contacts = output["contacts"][0, : len(crop), : len(crop)].float().cpu().numpy()
                inference_seconds = time.perf_counter() - inference_started
            else:
                for pair in pairs:
                    results.append(blank_result(locus, pair, len(sequence), args.model, str(exc)))
                continue

        for pair in pairs:
            hit_a = pair["hit_a"]
            hit_b = pair["hit_b"]
            a_start = int(hit_a["_start"]) - crop_start + 1
            a_end = int(hit_a["_end"]) - crop_start + 1
            b_start = int(hit_b["_start"]) - crop_start + 1
            b_end = int(hit_b["_end"]) - crop_start + 1
            row = blank_result(locus, pair, len(sequence), args.model, "")
            row.update(
                {
                    "crop_start": crop_start,
                    "crop_end": crop_end,
                    "crop_length": len(crop),
                    "device": str(device),
                    "inference_seconds": f"{inference_seconds:.4f}",
                }
            )
            try:
                metrics = contact_metrics(np, contacts, a_start, a_end, b_start, b_end)
                status, reason = classify_contact(metrics)
                row.update(
                    {
                        "max_contact_probability": f"{metrics['max']:.6f}",
                        "mean_contact_probability": f"{metrics['mean']:.6f}",
                        "top3_mean_contact_probability": f"{metrics['top3']:.6f}",
                        "contact_pairs_ge_0_5": metrics["n_ge_0_5"],
                        "separation_matched_percentile": f"{metrics['percentile']:.4f}",
                        "sequence_separation": metrics["separation"],
                        "contact_status": status,
                        "contact_reason": reason,
                    }
                )
            except Exception as exc:
                row["contact_status"] = "contact_unavailable"
                row["contact_reason"] = str(exc)
            results.append(row)

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    counts = defaultdict(int)
    for row in results:
        counts[row.get("contact_status", "unknown")] += 1
    summary = ", ".join(f"{key}={counts[key]}" for key in sorted(counts)) or "no applicable motif pairs"
    print(f"Wrote {output.resolve()} ({len(results)} pairs; {summary})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
