from __future__ import annotations

import json

import numpy as np

from witnesscell.cli import main

from ._fixture import synthetic_problem


def test_cli_fit_inspect_predict_validate(tmp_path, capsys) -> None:
    moments, split, gene2go, pairs = synthetic_problem()
    moments_path = tmp_path / "moments.npz"
    split_path = tmp_path / "split.json"
    go_path = tmp_path / "gene2go.json"
    model_path = tmp_path / "model.wcell"
    output_path = tmp_path / "prediction.npz"
    moments.to_npz(moments_path)
    split_path.write_text(
        json.dumps(
            {
                "train_conditions": split.train_conditions,
                "validation_conditions": split.validation_conditions,
            }
        ),
        encoding="utf-8",
    )
    go_path.write_text(json.dumps(gene2go), encoding="utf-8")

    assert main(["validate", "--moments", str(moments_path)]) == 0
    assert main(
        [
            "fit", "--moments", str(moments_path), "--split", str(split_path),
            "--gene2go", str(go_path), "--output", str(model_path),
        ]
    ) == 0
    assert main(["inspect", "--model", str(model_path)]) == 0
    assert main(
        [
            "predict", "--model", str(model_path), "--conditions", pairs[10], pairs[11],
            "--output", str(output_path),
        ]
    ) == 0
    with np.load(output_path, allow_pickle=False) as prediction:
        assert prediction["means"].shape == (2, len(moments.genes))
        assert set(prediction.files) == {
            "conditions", "genes", "means", "effects", "factorized_means", "factorized_effects"
        }
    assert '"valid": true' in capsys.readouterr().out
