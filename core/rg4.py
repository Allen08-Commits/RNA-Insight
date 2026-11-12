"""RNA G-quadruplex candidate screening utilities."""

from __future__ import annotations

import re

import pandas as pd


CANONICAL_RG4_PATTERN = re.compile(
    r"G{3,}[ACGU]{1,7}"
    r"G{3,}[ACGU]{1,7}"
    r"G{3,}[ACGU]{1,7}"
    r"G{3,}"
)


def composition_window_score(sequence: str) -> float:
    """
    Calculate a simple G/C composition-window score.

    G contributes +1, C contributes -1, and A/U contribute 0.
    This is a transparent heuristic, not a full G4Hunter reproduction.
    """
    sequence = sequence.upper().replace("T", "U")

    if not sequence:
        return 0.0

    return sum(
        1 if base == "G" else -1 if base == "C" else 0
        for base in sequence
    ) / len(sequence)


def scan_rg4(
    records: pd.DataFrame,
    window_size: int = 25,
    score_threshold: float = 0.4,
) -> pd.DataFrame:
    """
    Screen FASTA records for canonical rG4 motifs and composition windows.

    Returns one row per candidate region.
    """
    if "sequence_id" not in records.columns or "sequence" not in records.columns:
        raise ValueError("Input requires sequence_id and sequence columns.")

    if window_size < 1:
        raise ValueError("Window size must be positive.")

    if not -1 <= score_threshold <= 1:
        raise ValueError("Score threshold must be between -1 and 1.")

    results = []

    for _, record in records.iterrows():
        sequence_id = str(record["sequence_id"])
        sequence = str(record["sequence"]).upper().replace("T", "U")

        # Canonical motif candidates.
        for match in CANONICAL_RG4_PATTERN.finditer(sequence):
            start, end = match.span()
            candidate = match.group()

            results.append({
                "sequence_id": sequence_id,
                "start": start,
                "end": end,
                "candidate_sequence": candidate,
                "candidate_type": "canonical_motif",
                "score": composition_window_score(candidate),
                "length": len(candidate),
            })

        # Composition windows are a separate heuristic screen.
        if len(sequence) >= window_size:
            for start in range(0, len(sequence) - window_size + 1):
                end = start + window_size
                candidate = sequence[start:end]
                score = composition_window_score(candidate)

                if score >= score_threshold:
                    results.append({
                        "sequence_id": sequence_id,
                        "start": start,
                        "end": end,
                        "candidate_sequence": candidate,
                        "candidate_type": "composition_window",
                        "score": score,
                        "length": len(candidate),
                    })

    columns = [
        "sequence_id",
        "start",
        "end",
        "candidate_sequence",
        "candidate_type",
        "score",
        "length",
    ]

    if not results:
        return pd.DataFrame(columns=columns)

    result = pd.DataFrame(results)
    result = result.drop_duplicates(
        subset=["sequence_id", "start", "end", "candidate_type"]
    )

    return result.sort_values(
        ["sequence_id", "start", "end"]
    ).reset_index(drop=True)
