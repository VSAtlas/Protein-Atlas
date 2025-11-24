import logging
import os

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
    "vina.call",
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
    raw_topics = (os.environ.get("LOG_TOPICS") or str(cfg.get("LOG_TOPICS", ""))).replace(",", " ")
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
            tag = message[1:message.find("]")].strip().lower()
            allowed = self.allowed
            if not allowed or "all" in allowed:
                if (tag in DEFAULT_MUTED_TAGS) and ("all" not in allowed) and (tag not in allowed):
                    return False
                return True
            return tag in allowed

        lowered = message.lower()
        if "altloc" in lowered:
            tag = "altloc"
            allowed = self.allowed or set()
            if (not allowed or "all" in allowed) and (tag in DEFAULT_MUTED_TAGS) and ("all" not in allowed) and (tag not in allowed):
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
