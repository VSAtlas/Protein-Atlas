import atexit
import logging
import os
import sys
from pathlib import Path
from typing import Dict

from input_and_export_functions import _to_bool

# LOG_TOPICS semantics:
#   - By default (LOG_TOPICS unset/empty), INFO/DEBUG messages whose tag or
#     pseudo-tag is in DEFAULT_MUTED_TAGS are suppressed to keep logs readable.
#     This includes noisy internal topics (vina.call, cfg.emit, emit.debug,
#     propka.choice, read_any, elemfix), prep topics (ligprep, bulksdf, adt,
#     choose-mol2, arom-rescue, router.debug), and the pseudo-topic altloc
#     (altLoc/AltLocs found messages).
#   - Warnings and errors are always shown regardless of tag.
#   - Setting LOG_TOPICS=all shows all tagged and altLoc messages.
#   - Setting LOG_TOPICS to a whitespace/comma-separated list of tags
#     (e.g. "elemfix vina.call ligprep altloc") explicitly enables only those
#     topics (plus untagged messages, unless restricted elsewhere).
DEFAULT_MUTED_TAGS: set[str] = {
    # Very noisy internal debug topics
    # "vina.call",
    "cfg.emit",
    "emit.debug",
    "propka.choice",
    "read_any",
    "elemfix",
    # Element-repair noise from activesite/elem-fix
    "element",
    "elements",
    # Now also mute these by default (opt-in via LOG_TOPICS)
    "ligprep",
    "bulksdf",
    "adt",
    "choose-mol2",
    "arom-rescue",
    "router.debug",
    "RDKit",
    # Pseudo-topic for altLoc-style messages (handled specially)
    "altloc",
}


def coerce_log_level(name: str | None, default: int) -> int:
    """Convert a string level name to a logging level, falling back to default."""
    mapping = {
        "CRITICAL": logging.CRITICAL,
        "ERROR": logging.ERROR,
        "WARN": logging.WARNING,
        "WARNING": logging.WARNING,
        "INFO": logging.INFO,
        "DEBUG": logging.DEBUG,
        "NOTSET": logging.NOTSET,
    }
    return mapping.get(str(name or "").strip().upper(), default)


def parse_log_topics(cfg: dict) -> set[str]:
    """
    Return a lowercase set of tokens from LOG_TOPICS (env or cfg),
    splitting on commas and whitespace.
    """
    raw_topics = (
        os.environ.get("LOG_TOPICS") or str(cfg.get("LOG_TOPICS", ""))
    ).replace(",", " ")
    return {token.strip().lower() for token in raw_topics.split() if token.strip()}


class TopicFilter(logging.Filter):
    def __init__(self, allowed: set[str]):
        super().__init__()
        self.allowed = allowed

    def filter(self, record: logging.LogRecord) -> bool:
        if record.levelno >= logging.WARNING:
            return True
        message = record.getMessage()
        if message.startswith("[") and ("]" in message):
            tag = message[1 : message.find("]")].strip().lower()
            allowed = self.allowed
            if not allowed or "all" in allowed:
                if (
                    (tag in DEFAULT_MUTED_TAGS)
                    and ("all" not in allowed)
                    and (tag not in allowed)
                ):
                    return False
                return True
            return tag in allowed

        lowered = message.lower()
        if "altloc" in lowered:
            tag = "altloc"
            allowed = self.allowed or set()
            if (
                (not allowed or "all" in allowed)
                and (tag in DEFAULT_MUTED_TAGS)
                and ("all" not in allowed)
                and (tag not in allowed)
            ):
                return False
            if allowed and ("all" not in allowed):
                return tag in allowed
            return True

        return ("untagged" in self.allowed) or (not self.allowed)


def build_topic_filter(cfg: dict) -> logging.Filter | None:
    """
    Return a TopicFilter instance when LOG_TOPICS is non-empty and does not include 'all'.
    Return None when no topic filtering should be applied.
    """
    topics = parse_log_topics(cfg)
    if not topics or ("all" in topics):
        return None
    return TopicFilter(topics)


class _Tee:
    def __init__(self, stream, file_path):
        self._stream = stream
        self._fh = open(file_path, "a", buffering=1, encoding="utf-8", errors="replace")

    def write(self, data):
        try:
            self._stream.write(data)
        except Exception:
            pass
        try:
            self._fh.write(data)
        except Exception:
            pass

    def flush(self):
        try:
            self._stream.flush()
        except Exception:
            pass
        try:
            self._fh.flush()
        except Exception:
            pass

    def isatty(self):
        try:
            return self._stream.isatty()
        except Exception:
            return False

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass


def _tee_stdio_to(log_path):
    tee_out, tee_err = _Tee(sys.stdout, log_path), _Tee(sys.stderr, log_path)
    sys.stdout, sys.stderr = tee_out, tee_err

    def _announce_and_close():
        try:
            sys.stdout.write(f"Log -> {log_path}\n")
            sys.stdout.flush()
        finally:
            try:
                tee_out.close()
                tee_err.close()
            except Exception:
                pass

    atexit.register(_announce_and_close)


def make_protein_logger(docked_dir: str, pdb_id: str, cfg: Dict) -> logging.Logger:
    """
    Create a logger writing to DOCKED_DIR/<pdb_id>/protein.log and also to console.
    Keeps logs per-protein and avoids duplicate handlers.
    """

    base = Path(docked_dir)
    # If caller already passed .../docked/<PDB>, don't append <PDB> again
    log_dir = base if base.name.upper() == pdb_id.upper() else (base / pdb_id)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_dir / "protein.log"

    logger = logging.getLogger(pdb_id)
    logger.setLevel(logging.DEBUG)

    # Reset handlers to avoid duplicates if re-used
    if logger.hasHandlers():
        logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # Overwrite per run (truncate), not append
    fh = logging.FileHandler(log_file_path, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)

    ch = logging.StreamHandler(stream=sys.stdout)  # stdout, not stderr
    env_override = os.environ.get("QUIET_CONSOLE_OVERRIDE", "").strip()
    quiet = (
        (env_override.lower() in {"1", "true", "yes"})
        if env_override
        else _to_bool(cfg.get("QUIET_CONSOLE", False))
    )
    ch.setLevel(logging.WARNING if quiet else logging.INFO)
    ch.setFormatter(formatter)

    # Optional level overrides
    fh.setLevel(
        coerce_log_level(
            os.environ.get("LOG_LEVEL_FILE") or cfg.get("LOG_LEVEL_FILE"), fh.level
        )
    )
    ch.setLevel(
        coerce_log_level(
            os.environ.get("LOG_LEVEL_CONSOLE") or cfg.get("LOG_LEVEL_CONSOLE"),
            ch.level,
        )
    )

    # NOTE: logging_topics.TopicFilter mutes noisy topics (vina.call, ligprep, altloc, etc.)
    #       by default; use LOG_TOPICS=all or a list (e.g. LOG_TOPICS=ligprep,altloc) to opt back in.
    topic_filter = build_topic_filter(cfg)
    if topic_filter is not None:
        fh.addFilter(topic_filter)
        ch.addFilter(topic_filter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.propagate = False
    return logger


def bootstrap_root_logging(cfg: Dict, run_log_path: str) -> logging.Logger:
    """Ensure the root logger emits INFO-level records to the tee'd console."""

    root = logging.getLogger()
    if getattr(root, "_atlas_bootstrapped", False):
        return root

    # Clear any pre-existing handlers (e.g., from logging.basicConfig in imported modules)
    for h in list(root.handlers):
        root.removeHandler(h)
        try:
            h.close()
        except Exception:
            pass

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    configured_level = (
        os.environ.get("LOG_LEVEL_CONSOLE") or cfg.get("LOG_LEVEL_CONSOLE") or "INFO"
    )
    level_value = coerce_log_level(configured_level, logging.INFO)
    if level_value < logging.INFO:
        level_value = logging.INFO
    stream_handler.setLevel(level_value)

    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(name)s - %(message)s"
    )
    stream_handler.setFormatter(formatter)

    # NOTE: logging_topics.TopicFilter mutes noisy topics (vina.call, ligprep, altloc, etc.)
    #       by default; use LOG_TOPICS=all or a list (e.g. LOG_TOPICS=ligprep,altloc) to opt back in.
    topic_filter = build_topic_filter(cfg)
    if topic_filter is not None:
        stream_handler.addFilter(topic_filter)

    root.setLevel(logging.DEBUG)
    root.addHandler(stream_handler)
    root._atlas_bootstrapped = True  # type: ignore[attr-defined]
    return root
