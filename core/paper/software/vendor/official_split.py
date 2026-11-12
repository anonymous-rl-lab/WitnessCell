#!/usr/bin/env python3
"""Dependency-light reproduction of the official GEARS simulation split.

The code intentionally preserves NumPy's legacy RNG and the operation order in
GEARS 0.1.0.  The resulting Wessels test sets were checked condition-for-
condition against the published scPerturBench result tables for seeds 1--3.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pickle

import numpy as np


def load_gears_supported_genes(asset_dir: str | Path) -> frozenset[str]:
    """Reproduce the perturbation filter applied by GEARS 0.1.0.

    scPerturBench creates the published split *after* ``PertData.load`` removes
    perturbations absent from GEARS' default GO graph.  Omitting this step
    silently changes the random gene universe and therefore every downstream
    train/validation/test assignment.
    """
    asset_dir = Path(asset_dir)
    with (asset_dir / "gene2go_all.pkl").open("rb") as stream:
        gene2go = pickle.load(stream)
    with (asset_dir / "essential_all_data_pert_genes.pkl").open("rb") as stream:
        essential = pickle.load(stream)
    return frozenset(map(str, set(gene2go).intersection(essential)))


def filter_gears_supported_conditions(
    conditions: list[str], supported_genes: frozenset[str]
) -> list[str]:
    return [
        value
        for value in conditions
        if value == "control"
        or all(gene in supported_genes for gene in value.split("+"))
    ]


def clean_condition(condition: str) -> str:
    return condition.replace("+ctrl", "").replace("ctrl+", "").strip()


def gears_condition(condition: str) -> str:
    if condition == "control":
        return "ctrl"
    if "+" in condition:
        return condition
    return f"{condition}+ctrl"


def genes_from_perts(perts: list[str]) -> np.ndarray:
    return np.unique([
        gene
        for perturbation in np.unique(perts)
        for gene in perturbation.split("+")
        if gene != "ctrl"
    ])


def perts_from_genes(
    genes: np.ndarray, perts: list[str], kind: str = "both"
) -> list[str]:
    singles = [p for p in perts if "ctrl" in p and p != "ctrl"]
    combos = [p for p in perts if "ctrl" not in p]
    candidates = singles if kind == "single" else combos if kind == "combo" else perts
    return [
        perturbation
        for perturbation in candidates
        if any(gene in perturbation.split("+") for gene in genes)
    ]


def simulation_split(
    perts: list[str],
    train_gene_fraction: float,
    combo_seen2_train_fraction: float,
    seed: int,
) -> tuple[list[str], list[str], dict[str, list[str]]]:
    genes = genes_from_perts(perts)
    np.random.seed(seed)
    train_genes = np.random.choice(
        genes, int(len(genes) * train_gene_fraction), replace=False
    )
    ood_genes = np.setdiff1d(genes, train_genes)

    train = perts_from_genes(train_genes, perts, "single")
    combos = perts_from_genes(train_genes, perts, "combo")
    combo_seen1 = [
        condition
        for condition in combos
        if sum(gene in train_genes for gene in condition.split("+")) == 1
    ]
    combos = np.setdiff1d(combos, combo_seen1)
    np.random.seed(seed)
    combo_train = np.random.choice(
        combos, int(len(combos) * combo_seen2_train_fraction), replace=False
    )
    combo_seen2 = np.setdiff1d(combos, combo_train).tolist()
    unseen_single = perts_from_genes(ood_genes, perts, "single")
    ood_combos = perts_from_genes(ood_genes, perts, "combo")
    combo_seen0 = [
        condition
        for condition in ood_combos
        if sum(gene in train_genes for gene in condition.split("+")) == 0
    ]
    train.extend(combo_train.tolist())
    test = combo_seen1 + combo_seen2 + unseen_single + combo_seen0
    assert (
        len(combo_seen1)
        + len(combo_seen0)
        + len(unseen_single)
        + len(train)
        + len(combo_seen2)
        == len(perts)
    )
    return train, test, {
        "combo_seen0": combo_seen0,
        "combo_seen1": combo_seen1,
        "combo_seen2": combo_seen2,
        "unseen_single": unseen_single,
    }


@dataclass(frozen=True)
class OfficialSplit:
    train: tuple[str, ...]
    validation: tuple[str, ...]
    test: tuple[str, ...]
    test_subgroup: dict[str, tuple[str, ...]]
    validation_subgroup: dict[str, tuple[str, ...]]


def make_official_split(conditions: list[str], seed: int) -> OfficialSplit:
    unique = [gears_condition(value) for value in conditions]
    unique = list(dict.fromkeys(value for value in unique if value != "ctrl"))
    first_train, test, test_subgroup = simulation_split(unique, 0.8, 0.75, seed)
    train, validation, validation_subgroup = simulation_split(
        first_train, 0.9, 0.9, seed
    )
    train.append("ctrl")
    convert = lambda values: tuple(clean_condition(value) if value != "ctrl" else "control" for value in values)
    return OfficialSplit(
        train=convert(train),
        validation=convert(validation),
        test=convert(test),
        test_subgroup={key: convert(value) for key, value in test_subgroup.items()},
        validation_subgroup={key: convert(value) for key, value in validation_subgroup.items()},
    )
