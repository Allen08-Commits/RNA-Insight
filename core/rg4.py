
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
    Assign run-based G4Hunter scores.

    Consecutive G bases receive positive scores.
    Consecutive C bases receive negative scores.
    Run scores are capped at 4.
    A and U receive zero.
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
    """Merge overlapping G4Hunter-positive sliding windows."""
    if not hits:
        return []

    hits.sort(key=lambda item: (item["start"], item["end"]))

    groups = []
    current = [hits[0]]

    for hit in hits[1:]:
        current_end = max(item["end"] for item in current)

        if hit["start"] < current_end:
            current.append(hit)
        else:
            groups.append(current)
            current = [hit]

    groups.append(current)

    merged = []

    for group in groups:
        start = min(item["start"] for item in group)
        end = max(item["end"] for item in group)

        candidate = sequence[start:end]
        candidate_scores = calculate_g4hunter_scores(candidate)

        merged.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": candidate,
            "candidate_type": "g4hunter",
            "score": round(
                sum(candidate_scores) / len(candidate_scores), 4
            ),
            "length": len(candidate),
            "n_windows": len(group),
            "method": "G4Hunter-style",
        })

    return merged


def _find_canonical_motifs(
    sequence: str,
    sequence_id: str,
) -> list[dict]:
    """Identify canonical four-tract G-rich RNA motif candidates."""
    results = []

    for match in CANONICAL_RG4_PATTERN.finditer(sequence):
        candidate = match.group(1)
        start = match.start(1)
        end = match.end(1)

        scores = calculate_g4hunter_scores(candidate)

        results.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": candidate,
            "candidate_type": "canonical_motif",
            "score": round(sum(scores) / len(scores), 4),
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
    Identify RNA G4 candidate regions using:

    1. G4Hunter-style sliding-window scoring.
    2. Independent canonical four-tract motif screening.

    Overlapping positive windows are merged into candidate regions.
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

        # G4Hunter-style sliding-window analysis.
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
