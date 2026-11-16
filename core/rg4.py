
"""RNA G-quadruplex prediction using G4Hunter-style scoring."""

from __future__ import annotations

import re

import pandas as pd


CANONICAL_RG4_PATTERN = re.compile(
    r"(?=(G{3,}[ACGU]{1,7}G{3,}[ACGU]{1,7}"
    r"G{3,}[ACGU]{1,7}G{3,}))"
)


def calculate_g4hunter_scores(sequence: str) -> list[int]:
    """
    Assign per-nucleotide G4Hunter scores.

    G runs receive positive scores and C runs negative scores.
    The magnitude is the run length, capped at 4.
    A and U receive zero.

    This function returns scores in transcript orientation.
    """
    sequence = sequence.upper().replace("T", "U")

    scores = [0] * len(sequence)
    index = 0

    while index < len(sequence):
        base = sequence[index]

        if base not in {"G", "C"}:
            index += 1
            continue

        end = index + 1

        while end < len(sequence) and sequence[end] == base:
            end += 1

        run_length = end - index
        magnitude = min(run_length, 4)

        score = magnitude if base == "G" else -magnitude

        for position in range(index, end):
            scores[position] = score

        index = end

    return scores


def _merge_window_hits(
    sequence: str,
    hits: list[dict],
    sequence_id: str,
) -> list[dict]:
    """Merge overlapping G4Hunter sliding-window hits."""
    if not hits:
        return []

    hits.sort(key=lambda hit: (hit["start"], hit["end"]))

    groups = []
    current = [hits[0]]

    for hit in hits[1:]:
        if hit["start"] < max(item["end"] for item in current):
            current.append(hit)
        else:
            groups.append(current)
            current = [hit]

    groups.append(current)

    merged = []

    for group in groups:
        start = min(item["start"] for item in group)
        end = max(item["end"] for item in group)

        region_scores = calculate_g4hunter_scores(
            sequence[start:end]
        )

        score = sum(region_scores) / len(region_scores)

        merged.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": sequence[start:end],
            "candidate_type": "g4hunter",
            "score": round(score, 4),
            "length": end - start,
            "n_windows": len(group),
            "method": "G4Hunter-style",
        })

    return merged


def _find_canonical_motifs(
    sequence: str,
    sequence_id: str,
) -> list[dict]:
    """Find canonical four-tract RNA G4 motif patterns."""
    results = []

    for match in CANONICAL_RG4_PATTERN.finditer(sequence):
        candidate = match.group(1)
        start = match.start(1)
        end = match.end(1)

        results.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": candidate,
            "candidate_type": "canonical_motif",
            "score": round(
                sum(calculate_g4hunter_scores(candidate))
                / len(candidate),
                4,
            ),
            "length": len(candidate),
            "n_windows": 0,
            "method": "Canonical G-tract motif",
        })

    return results


def scan_rg4(
    records: pd.DataFrame,
    window_size: int = 25,
    score_threshold: float = 1.2,
) -> pd.DataFrame:
    """
    Identify RNA G4 candidates using G4Hunter-style scoring
    and an independent canonical motif screen.

    Only positive G-rich windows are reported as RNA candidates.
    Overlapping qualifying windows are merged into regions.
    """
    required = {"sequence_id", "sequence"}

    if not required.issubset(records.columns):
        raise ValueError(
            "Input requires sequence_id and sequence columns."
        )

    if window_size < 1:
        raise ValueError("Window size must be positive.")

    if not 0 <= score_threshold <= 4:
        raise ValueError("Score threshold must be between 0 and 4.")

    output = []

    for _, record in records.iterrows():
        sequence_id = str(record["sequence_id"])
        sequence = str(record["sequence"]).upper().replace("T", "U")

        if not sequence:
            continue

        scores = calculate_g4hunter_scores(sequence)

        # G4Hunter sliding-window analysis.
        if len(sequence) >= window_size:
            hits = []

            for start in range(len(sequence) - window_size + 1):
                end = start + window_size
                window_score = (
                    sum(scores[start:end]) / window_size
                )

                if window_score >= score_threshold:
                    hits.append({
                        "start": start,
                        "end": end,
                        "window_score": window_score,
                    })

            output.extend(
                _merge_window_hits(
                    sequence,
                    hits,
                    sequence_id,
                )
            )

        # Independent canonical motif analysis.
        output.extend(
            _find_canonical_motifs(sequence, sequence_id)
        )

    columns = [
        "sequence_id",
        "start",
        "end",
        "candidate_sequence",
        "candidate_type",
        "score",
        "length",
        "n_windows",
        "method",
    ]

    if not output:
        return pd.DataFrame(columns=columns)

    result = pd.DataFrame(output, columns=columns)

    return result.sort_values(
        ["sequence_id", "start", "end", "candidate_type"]
    ).reset_index(drop=True)
