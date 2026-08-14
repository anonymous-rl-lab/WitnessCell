from __future__ import annotations

import argparse
import ast
import importlib.util
import os
import sys
import types
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.metric_core import (  # noqa: E402
    drf,
    nir,
    source_weight_transform,
    weighted_delta_r2,
    wmse,
)


def _load_upstream(normative_repo: Path):
    package = types.ModuleType("e22_locked_cellsimbench")
    package.__path__ = []
    core = types.ModuleType("e22_locked_cellsimbench.core")
    core.__path__ = []
    sys.modules[package.__name__] = package
    sys.modules[core.__name__] = core

    def load(name: str, path: Path):
        spec = importlib.util.spec_from_file_location(name, path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot load locked source {path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    data_manager = load(
        "e22_locked_cellsimbench.core.data_manager",
        normative_repo / "src/cellsimbench/core/data_manager.py",
    )
    metrics_engine = load(
        "e22_locked_cellsimbench.core.metrics_engine",
        normative_repo / "src/cellsimbench/core/metrics_engine.py",
    )
    return data_manager, metrics_engine


def _extract_locked_drf(normative_repo: Path):
    path = normative_repo / "analyses/calibration/calibration_analysis.py"
    tree = ast.parse(path.read_text())
    candidates = [node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name == "calculate_drf"]
    if len(candidates) != 1:
        raise RuntimeError(f"expected one nested calculate_drf, found {len(candidates)}")
    function = candidates[0]
    function.decorator_list = []
    module = ast.Module(body=[ast.Import(names=[ast.alias(name="numpy", asname="np")]), function], type_ignores=[])
    ast.fix_missing_locations(module)
    namespace: dict[str, object] = {}
    exec(compile(module, str(path), "exec"), namespace)
    return namespace["calculate_drf"]


class SourceParityTests(unittest.TestCase):
    normative_repo: Path

    @classmethod
    def setUpClass(cls) -> None:
        value = os.environ.get("E22_NORMATIVE_REPO")
        if not value:
            raise RuntimeError("E22_NORMATIVE_REPO is required for source-parity tests")
        cls.normative_repo = Path(value)
        cls.up_data, cls.up_metrics = _load_upstream(cls.normative_repo)

    def test_wmse_and_weighted_r2(self) -> None:
        pred = np.array([0.2, -1.0, 2.5, 0.0])
        truth = np.array([0.0, -0.5, 2.0, 1.0])
        baseline = np.array([0.1, 0.1, 0.1, 0.1])
        weight = source_weight_transform([1, 4, 2, 3], ["a", "b", "c", "d"], ["a", "b", "c", "d"])
        self.assertEqual(wmse(pred, truth, weight), self.up_data.wmse(pred, truth, weight.values))
        self.assertEqual(
            weighted_delta_r2(pred, truth, baseline, weight),
            self.up_data.r2_score_on_deltas(truth - baseline, pred - baseline, weight.values),
        )

    def test_weight_transform(self) -> None:
        manager = self.up_data.DataManager.__new__(self.up_data.DataManager)
        manager.deg_scores_dict = {"x_target": np.array([-1.0, 3.0, 2.0, 5.0])}
        manager.deg_names_dict = {"x_target": np.array(["g1", "g2", "g2", "unused"])}
        manager.pert_normalized_abs_scores_vsrest = {}
        manager.adata = type("Adata", (), {"var_names": pd.Index(["g2", "g1", "g3"])})()
        manager._precompute_deg_weights()
        local = source_weight_transform(
            manager.deg_scores_dict["x_target"], manager.deg_names_dict["x_target"], manager.adata.var_names
        )
        np.testing.assert_array_equal(local.values, manager.pert_normalized_abs_scores_vsrest["x_target"])

    def test_nir(self) -> None:
        identities = ["cov_a", "cov_b", "cov_c"]
        predictions = np.array([[0.1, 0.0], [1.0, 0.2], [2.0, -0.1]])
        truths = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
        frame_pred = pd.DataFrame(predictions, index=identities)
        frame_true = pd.DataFrame(truths, index=identities)
        engine = self.up_metrics.MetricsEngine.__new__(self.up_metrics.MetricsEngine)
        upstream = engine._calculate_nir_scores(frame_pred, frame_true)
        self.assertEqual(nir(predictions, truths, identities), upstream)

    def test_drf_ast_extracted_source(self) -> None:
        upstream = _extract_locked_drf(self.normative_repo)
        for baseline, duplicate, higher in [(4.0, 1.0, False), (0.2, 0.8, True), (1.0, 3.0, False)]:
            config = {"higher_better": higher, "perfect": 1.0 if higher else 0.0}
            self.assertEqual(drf(baseline, duplicate, higher_better=higher), upstream(baseline, duplicate, config))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--normative-repo", type=Path, required=True)
    args, remaining = parser.parse_known_args(sys.argv[1:])
    os.environ["E22_NORMATIVE_REPO"] = str(args.normative_repo)
    unittest.main(argv=[sys.argv[0], *remaining])
