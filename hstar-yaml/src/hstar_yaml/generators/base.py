"""Base class for HSTAR file generators."""

from __future__ import annotations

from pathlib import Path

from ..schema import HstarProject
from ..mesh_reader import MeshData


class BaseGenerator:
    """Base class for generating HSTAR input files."""

    extension: str = ""

    def __init__(self, project: HstarProject, mesh: MeshData | None = None):
        self.project = project
        self.mesh = mesh

    def generate(self, output_dir: Path) -> Path | None:
        """Generate the file. Returns the path or None if not applicable."""
        content = self.build()
        if content is None:
            return None
        path = output_dir / f"{self.project.problem_name}{self.extension}"
        path.write_text(content, encoding="utf-8")
        return path

    def build(self) -> str | None:
        """Build file content as string. Override in subclasses."""
        raise NotImplementedError

    def _fmt(self, value: float, width: int = 14) -> str:
        """Format a float value for Fortran free-format read."""
        # Use scientific notation for very large/small values
        if value == 0.0:
            return f"{0.0:{width}.3E}"
        abs_val = abs(value)
        if abs_val >= 1e7 or abs_val < 1e-2:
            return f"{value:{width}.3E}"
        return f"{value:{width}.6f}"

    def _fmtg(self, value: float) -> str:
        """Format a float in Fortran-style compact format.

        Rules matching HSTAR convention:
        - Very small/large values: use E notation like 0.240E+04
        - Integers: use plain int format
        - Normal range: use compact decimal
        """
        if value == 0.0:
            return "0.0"
        abs_val = abs(value)
        if abs_val >= 1e7 or abs_val < 1e-3:
            return f"{value:.3E}"
        if value == int(value) and abs_val < 1e6:
            return f"{int(value)}.0"
        return f"{value:g}"
