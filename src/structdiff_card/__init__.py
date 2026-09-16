"""StructDiff Card public API."""

from .analysis import compare_structures
from .models import (
    AlignmentResult,
    ChainAlignment,
    ContactChange,
    ContactSummary,
    DifferenceRegion,
    DifferenceSummary,
    MotionDecomposition,
    ResidueMapping,
    SpatialPatch,
    Transform,
)

__all__ = [
    "AlignmentResult",
    "ChainAlignment",
    "ContactChange",
    "ContactSummary",
    "DifferenceRegion",
    "DifferenceSummary",
    "MotionDecomposition",
    "ResidueMapping",
    "SpatialPatch",
    "Transform",
    "compare_structures",
]

__version__ = "0.4.0"
