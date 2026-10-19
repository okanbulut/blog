"""Read and prepare the PERSUADE 2.0 essays for automated essay scoring.

Companion code for the blog post "DeBERTa vs. UltraBERT: Does a Longer Context Window
Help Automated Essay Scoring?" by Okan Bulut.

Before running:
1. Download the two data files from the links in the PERSUADE 2.0 README:
   https://github.com/scrosseye/persuade_corpus_2.0
   The test set is a password-protected zip file; the password is in the README.
   The data are licensed under CC BY-NC-SA 4.0.
2. Put persuade_corpus_2.0_train.csv and persuade_corpus_2.0_test.csv in the same folder
   as this script (or change DATA_DIR below).

Run:
    python prepare_persuade.py

Output:
    persuade_essays.csv  one row per essay: essay_id, text, score, prompt_name, task,
                         ell_status, split (train / dev / test)

Requires: pandas, numpy
"""

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(".")
FILES = ["persuade_corpus_2.0_train.csv", "persuade_corpus_2.0_test.csv"]
SEED = 42


def read_essays(data_dir: Path) -> pd.DataFrame:
    """Read both files and keep one row per essay.

    The files have one row per *discourse element* (lead, position, claim, ...), and the
    essay-level columns are repeated on every row of the same essay.
    """
    rows = pd.concat([pd.read_csv(data_dir / f, low_memory=False) for f in FILES],
                     ignore_index=True)
    print("Rows (discourse elements):", rows.shape)
    essays = (rows[["essay_id_comp", "full_text", "holistic_essay_score",
                    "prompt_name", "task", "ell_status"]]
              .drop_duplicates("essay_id_comp")
              .rename(columns={"essay_id_comp": "essay_id",
                               "full_text": "text",
                               "holistic_essay_score": "score"})
              .reset_index(drop=True))
    print("Essays:", essays.shape)
    return essays


def text_id(text: str) -> str:
    """A stable ID made from the essay text itself ("X" + 11 hex characters)."""
    key = re.sub(r"\s+", " ", text).strip().lower()
    return "X" + hashlib.sha1(key.encode()).hexdigest()[:11].upper()


def clean_essays(essays: pd.DataFrame) -> pd.DataFrame:
    essays = essays.copy()

    # Problem 1: IDs damaged by spreadsheet software (e.g., "4.45649E+11").
    # Real IDs are 12 uppercase hexadecimal characters.
    valid = essays["essay_id"].str.fullmatch(r"[0-9A-F]{12}")
    print("Damaged IDs:", (~valid).sum())
    essays.loc[~valid, "essay_id"] = essays.loc[~valid, "text"].map(text_id)

    # Problem 2: the same essay text appearing under two different IDs.
    key = essays["text"].str.replace(r"\s+", " ", regex=True).str.strip().str.lower()
    print("Duplicate texts:", key.duplicated().sum())
    essays = essays[~key.duplicated()].reset_index(drop=True)

    # Problem 3: blank strings that actually mean "missing".
    essays["ell_status"] = essays["ell_status"].str.strip().replace("", np.nan)
    return essays


def stratified_split(df: pd.DataFrame, dev_frac: float = 0.1, test_frac: float = 0.2,
                     seed: int = SEED) -> pd.Series:
    """Shuffle each prompt x score cell and cut it into test, dev, and train."""
    rng = np.random.default_rng(seed)
    split = pd.Series("train", index=df.index)
    for _, cell in df.groupby(["prompt_name", "score"]):
        idx = rng.permutation(cell.index.to_numpy())
        n_test = round(len(idx) * test_frac)
        n_dev = round(len(idx) * dev_frac)
        split[idx[:n_test]] = "test"
        split[idx[n_test:n_test + n_dev]] = "dev"
    return split


if __name__ == "__main__":
    essays = clean_essays(read_essays(DATA_DIR))
    print("Essays after cleaning:", essays.shape)
    print(essays["score"].value_counts().sort_index())

    essays["split"] = stratified_split(essays)
    print(essays["split"].value_counts())

    essays.to_csv(DATA_DIR / "persuade_essays.csv", index=False)
    print("Saved persuade_essays.csv")
