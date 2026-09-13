"""Shared runtime configuration for PathwayLGM entry points."""

import logging
from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class RuntimeConfig:
    """Runtime settings shared by package commands.

    Parameters
    ----------
    log_level : str, default="INFO"
        Minimum severity for the package logger. Model-specific configuration
        will be introduced alongside each analysis step.
    """

    log_level: Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"] = "INFO"

    def __post_init__(self) -> None:
        """Reject unsupported logging levels before configuring handlers."""
        if self.log_level not in {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}:
            raise ValueError(f"Unsupported log level: {self.log_level}")


def configure_logging(config: RuntimeConfig) -> logging.Logger:
    """Configure the package logger without changing the root logger.

    Repeated calls reuse existing handlers. Imports alone do not configure
    logging; command-line entry points opt in explicitly.
    """
    logger = logging.getLogger("pathwaylgm")
    logger.setLevel(config.log_level)
    logger.propagate = False

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
        )
        logger.addHandler(handler)
    return logger
