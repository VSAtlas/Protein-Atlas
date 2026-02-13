from __future__ import annotations

import numpy as np
from scipy import sparse
from sklearn.linear_model import LogisticRegression

from ml.config import ModelConfig


def train_logistic_regression(
    *,
    X_train: sparse.csr_matrix,
    y_train: np.ndarray,
    model_config: ModelConfig,
    random_seed: int,
) -> LogisticRegression:
    if model_config.type.strip().lower() != "logreg":
        raise ValueError(f"Unsupported model.type={model_config.type!r}. Expected 'logreg'.")

    unique_labels = np.unique(y_train)
    if unique_labels.size < 2:
        raise ValueError(
            f"Logistic regression requires both classes. Found labels={unique_labels.tolist()}."
        )

    model = LogisticRegression(
        C=float(model_config.C),
        class_weight=model_config.class_weight,
        max_iter=int(model_config.max_iter),
        random_state=int(random_seed),
        solver="liblinear",
    )
    model.fit(X_train, y_train)
    return model

