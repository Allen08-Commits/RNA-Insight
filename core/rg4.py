
"""RNA G-quadruplex prediction and candidate classification.

Includes:
1. G4Hunter-style sliding-window screening.
2. Canonical four-tract G-rich motif screening.
3. Separate score fields for each prediction method.
4. Merging of overlapping G4Hunter-positive windows.
"""

from __future__ import annotations

import re

import pandas as pd


CANONICAL_RG4_PATTERN = re.compile(
    r"(?=(G{3,}[ACGU]{1,7}G{3,}[ACGU]{1,7}"
    r"G{3,}[ACGU]{1,7}G{3,}))"
)


RESULT_COLUMNS = [
    "sequence_id",
    "start",
    "end",
    "candidate_sequence",
    "candidate_type",
    "score",
    "score_type",
    "length",
    "n_windows",
    "method",
]


def calculate_g4hunter_scores(sequence: str) -> list[int]:
    """Calculate run-based G4Hunter-style nucleotide scores."""
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

        value = magnitude if base == "G" else -magnitude

        for position in range(index, end):
            scores[position] = value

        index = end

    return scores


def calculate_window_score(scores: list[int]) -> float:
    """Calculate the mean nucleotide score for a window."""
    if not scores:
        return 0.0

    return sum(scores) / len(scores)


def _merge_window_hits(
    sequence: str,
    hits: list[dict],
    sequence_id: str,
) -> list[dict]:
    """Merge overlapping G4Hunter-positive windows."""
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

    results = []

    for group in groups:
        start = min(item["start"] for item in group)
        end = max(item["end"] for item in group)

        candidate = sequence[start:end]

        # Average of the qualifying sliding-window scores.
        mean_window_score = sum(
            item["window_score"] for item in group
        ) / len(group)

        results.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": candidate,
            "candidate_type": "g4hunter",
            "score": round(mean_window_score, 4),
            "score_type": "Mean G4Hunter window score",
            "length": end - start,
            "n_windows": len(group),
            "method": "G4Hunter-style sliding window",
        })

    return results


def _find_canonical_motifs(
    sequence: str,
    sequence_id: str,
) -> list[dict]:
    """Find canonical four-tract G-rich RNA motif candidates."""
    results = []

    for match in CANONICAL_RG4_PATTERN.finditer(sequence):
        candidate = match.group(1)
        start = match.start(1)
        end = match.end(1)

        scores = calculate_g4hunter_scores(candidate)

        # This is a descriptive composition score, NOT a
        # G4Hunter sliding-window score.
        composition_score = calculate_window_score(scores)

        results.append({
            "sequence_id": sequence_id,
            "start": start,
            "end": end,
            "candidate_sequence": candidate,
            "candidate_type": "canonical_motif",
            "score": round(composition_score, 4),
            "score_type": "Motif composition score",
            "length": len(candidate),
            "n_windows": 0,
            "method": "Canonical four-tract G-rich motif",
        })

    return results


def scan_rg4(
    records: pd.DataFrame,
    window_size: int = 25,
    score_threshold: float = 1.2,
) -> pd.DataFrame:
    """
    Screen RNA sequences using two independent approaches.

    G4Hunter-style:
        Positive G-rich sliding-window scores.

    Canonical motif:
        Four G-rich tracts separated by short loop regions.

    The two methods are reported separately. Their intervals may overlap
    and should not be interpreted as independent structures.
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

        if set(sequence) - set("ACGUNRYKMSWBDHV"):
            continue

        nucleotide_scores = calculate_g4hunter_scores(sequence)

        # G4Hunter-style sliding-window analysis.
        if len(sequence) >= window_size:
            hits = []

            for start in range(len(sequence) - window_size + 1):
                end = start + window_size

                window_score = calculate_window_score(
                    nucleotide_scores[start:end]
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

    if not output:
        return pd.DataFrame(columns=RESULT_COLUMNS)

    result = pd.DataFrame(output, columns=RESULT_COLUMNS)

    return result.sort_values(
        ["sequence_id", "start", "end", "candidate_type"]
    ).reset_index(drop=True)
