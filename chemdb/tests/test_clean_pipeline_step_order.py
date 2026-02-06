
import protein_prep.clean_pipeline as clean_pipeline


def test_clean_pipeline_step_order_smoke(monkeypatch, tmp_path):
    call_order = []
    expected = [
        "init_workspace_and_paths",
        "step_altloc_filter",
        "step_strip_nonstandard",
        "step_modeller_loopfill",
        "step_phenix_clean",
        "step_reduce_protonation",
        "step_prepare_receptor",
        "finalize_outputs_and_logs",
    ]

    def _mk_step(name):
        def _step(ctx):
            call_order.append(name)
            return ctx

        return _step

    out_file = tmp_path / "cleaned.pdb"

    monkeypatch.setattr(
        clean_pipeline, "init_workspace_and_paths", _mk_step("init_workspace_and_paths")
    )
    monkeypatch.setattr(
        clean_pipeline, "step_altloc_filter", _mk_step("step_altloc_filter")
    )
    monkeypatch.setattr(
        clean_pipeline, "step_strip_nonstandard", _mk_step("step_strip_nonstandard")
    )
    monkeypatch.setattr(
        clean_pipeline, "step_modeller_loopfill", _mk_step("step_modeller_loopfill")
    )
    monkeypatch.setattr(
        clean_pipeline, "step_phenix_clean", _mk_step("step_phenix_clean")
    )
    monkeypatch.setattr(
        clean_pipeline, "step_reduce_protonation", _mk_step("step_reduce_protonation")
    )

    def _prepare(ctx):
        call_order.append("step_prepare_receptor")
        out_file.write_text("ATOM\n", encoding="utf-8")
        ctx["result"] = str(out_file)
        return ctx

    monkeypatch.setattr(clean_pipeline, "step_prepare_receptor", _prepare)

    def _finalize(ctx):
        call_order.append("finalize_outputs_and_logs")
        return ctx.get("result")

    monkeypatch.setattr(clean_pipeline, "finalize_outputs_and_logs", _finalize)

    result = clean_pipeline.clean_pdb("input.pdb", str(tmp_path))

    assert result == str(out_file)
    assert out_file.exists()
    assert call_order == expected
