from __future__ import annotations

from typing import Any


def heatmap_call(attr: str, /, *args: Any, **kwargs: Any) -> Any:
    import analysis.reporting.heatmap_html as surface

    return getattr(surface, attr)(*args, **kwargs)
