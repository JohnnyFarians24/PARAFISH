"""Central logging setup utilities."""

from __future__ import annotations

import logging


def configure_logging(name: str = "parasite_platform") -> logging.Logger:
	"""Create and return a project logger with a basic stream handler."""
	logger = logging.getLogger(name)
	if not logger.handlers:
		handler = logging.StreamHandler()
		formatter = logging.Formatter(
			"%(asctime)s | %(name)s | %(levelname)s | %(message)s"
		)
		handler.setFormatter(formatter)
		logger.addHandler(handler)
	logger.setLevel(logging.INFO)
	return logger
