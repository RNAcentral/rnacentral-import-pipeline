# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import polars as pl
import onnxruntime as rt
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)

def evaluate_model(model_path, test_data, exclude) -> pl.DataFrame:
    sess = rt.InferenceSession(model_path, providers=["CPUExecutionProvider"])
    input_name = sess.get_inputs()[0].name
    label_name = sess.get_outputs()[0].name
    probability_name = sess.get_outputs()[1].name

    test = pl.read_parquet(test_data)
    excludes = list(exclude)
    excludes.append("label")

    X_test = test.select(pl.exclude(excludes)).to_numpy()
    y_test = test.select("label").to_numpy().flatten()

    if X_test.height > 10_000:
        predictions = []
        prob_dict = []
        chunk_size = 1000
        for batch_idx in range(0, X_test.height, chunk_size):
            pred_part, prob_part = sess.run(
                [label_name, probability_name],
                {input_name: X_test[batch_idx : batch_idx + chunk_size]},
            )
            predictions.extend(pred_part)
            prob_dict.extend(prob_part)
    else:
        predictions, prob_dict = sess.run(
            [label_name, probability_name], {input_name: X_test}
        )

    probabilities = [pr[c] for pr, c in zip(prob_dict, predictions)]

    bal_acc = balanced_accuracy_score(y_test, predictions)
    f1 = f1_score(y_test, predictions)
    auc = roc_auc_score(y_test, probabilities)
    ap = average_precision_score(y_test, probabilities)

    metrics = {"balanced_acc": bal_acc, "F1": f1, "auc": auc, "ap": ap}
    metrics = pl.DataFrame(metrics)

    return metrics
