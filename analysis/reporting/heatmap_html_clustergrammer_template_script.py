from __future__ import annotations

from analysis.reporting.heatmap_html_clustergrammer_template_script_bootstrap import (
    build_clustergrammer_client_script_bootstrap,
)
from analysis.reporting.heatmap_html_clustergrammer_template_script_interactions import (
    build_clustergrammer_client_script_interactions,
)
from analysis.reporting.heatmap_html_clustergrammer_template_script_runtime import (
    build_clustergrammer_client_script_runtime,
)


def build_clustergrammer_client_script(**kwargs) -> str:
    return (
        build_clustergrammer_client_script_bootstrap(**kwargs)
        + build_clustergrammer_client_script_interactions(**kwargs)
        + build_clustergrammer_client_script_runtime(**kwargs)
    )
