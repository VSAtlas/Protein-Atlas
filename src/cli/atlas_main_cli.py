from __future__ import annotations


def main() -> int:
    # Compatibility wrapper only.
    # Keep user-facing QoL routing in cli.qol_cli and runtime behavior in cli.main_pipeline.
    # Defer import so CLI startup does not force heavy module import unless invoked.
    from cli.qol_cli import dispatch

    routed = dispatch()
    if routed is not None:
        return int(routed)

    from main import main as root_main

    root_main()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
