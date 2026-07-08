from __future__ import annotations

import random
from typing import Sequence, TypeVar

T = TypeVar("T")


def shuffled(values: Sequence[T], *, seed: int) -> list[T]:
    out = list(values)
    random.Random(seed).shuffle(out)
    return out

