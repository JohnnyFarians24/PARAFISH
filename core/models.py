"""Core contracts shared by every stage in the PARAFISH pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional

from core.config import MIN_UNIQUE_SIZE


class AnalysisMode(Enum):
	"""Library layout for input FASTQ files."""

	SE = "SE"
	PE = "PE"


class AnalysisBranch(Enum):
	"""Separated analysis branches."""

	VSEARCH = "vsearch"
	DADA2 = "dada2"


class PhyloMode(Enum):
	"""Phylogenetic routing modes."""

	FAST = "fasttree"
	ROBUST = "raxml"
	HYBRID = "hybrid"


class Marker(Enum):
	"""Marker families used for downstream routing and references."""

	COI = "COI"
	ITS = "ITS"
	SSU_18S = "18S"
	LSU_28S = "28S"
	OTHER = "OTHER"


@dataclass
class SamplePaths:
	"""Caminhos dos ficheiros FASTQ de uma amostra."""

	r1: Path
	r2: Optional[Path] = None

	def validate(self, mode: AnalysisMode):
		"""Validate that the requested layout exists on disk."""
		assert self.r1.exists(), f"R1 não encontrado: {self.r1}"
		if mode == AnalysisMode.PE:
			assert self.r2 is not None, "Modo PE requer R2"
			assert self.r2.exists(), f"R2 não encontrado: {self.r2}"


@dataclass
class Sample:
	"""Representa uma amostra individual."""

	internal_id: str
	display_label: str
	paths: SamplePaths
	mode: AnalysisMode
	marker: Marker
	selected: bool = True
	# Negative/blank controls are sequenced alongside real samples to expose reagent
	# and cross-contamination. They must never be pooled into diversity statistics;
	# they are the reference against which contaminant features are called.
	is_control: bool = False
	# Free-form experimental group (site, host species, treatment...). Drives the
	# grouped comparisons (PERMANOVA, ordination ellipses). None = ungrouped.
	group: Optional[str] = None


@dataclass
class StepOutput:
	"""Output padronizado de uma etapa do pipeline."""

	step_name: str
	success: bool
	files: Dict[str, Path]
	log: Optional[Path] = None
	metrics: Dict[str, float] = field(default_factory=dict)
	error_message: Optional[str] = None


@dataclass
class RunConfig:
	"""Configuração completa de um run."""

	run_id: str
	samples: List[Sample]
	branch: AnalysisBranch
	analysis_method: str = "vsearch"
	enabled_stages: Dict = field(default_factory=dict)
	phylo_mode: Optional[PhyloMode] = None
	output_dir: Path = Path("./output")
	pear_params: Dict = field(default_factory=lambda: {
		"min_overlap": 10,
		"min_length": 200,
		"max_length": 500,
		"p_value": 0.01,
		"max_uncalled_base": 0.0,
		"threads": 4,
		"test": 1,
	})
	trimmomatic_params: Dict = field(default_factory=lambda: {
		"leading": 3,
		"trailing": 3,
		"window_size": 4,
		"window_quality": 20,
		"minlen": 200,
	})
	vsearch_params: Dict = field(default_factory=lambda: {
		"fastq_maxee": 1.0,
		"fastq_minlen": 200,
		"id": 0.97,
		"sintax_cutoff": 0.8,
		# The smallest cluster allowed to SEED an OTU. Configurable via
		# PARAFISH_MIN_UNIQUE_SIZE; see core.config.MIN_UNIQUE_SIZE for the default (1)
		# and the measurements behind it. Not a literal here: the value is a scientific
		# choice a study makes, not a constant of the code.
		"min_unique_size": MIN_UNIQUE_SIZE,
		"threads": 8,
	})
	dada2_params: Dict = field(default_factory=lambda: {
		"truncLen_F": 0,
		"truncLen_R": 0,
		"maxEE_F": 2,
		"maxEE_R": 2,
		"truncQ": 2,
		"ref_fasta": None,
		"multithread": True,
	})
	mafft_params: Dict = field(default_factory=lambda: {
		"mode": "auto",
		"adjust_direction": True,
		"threads": 8,
	})
	phylo_params: Dict = field(default_factory=lambda: {
		"model": "GTR+G",
		"bs_metric": "fbp",
		"bs_trees": 10,
		"threads": 8,
		"seed": 42,
		"fasttree_gamma": True,
		"fasttree_gtr": True,
	})
	bold_params: Dict = field(default_factory=lambda: {
		"db": 1,
		"mode": 2,
		"thresholds": [97, 95, 90, 85],
	})


@dataclass
class DatabaseConfig:
	"""Database selection settings derived from marker and source."""

	marker: str
	source: str
	reference_path: Optional[Path] = None

