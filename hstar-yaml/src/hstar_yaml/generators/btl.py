"""Generate the .btl (back-analysis control) file.

Format:
    Line 1: nparams
    Line 2+: param definitions (material, property, initial, bounds)
"""

from __future__ import annotations

from .base import BaseGenerator


class BtlGenerator(BaseGenerator):
    extension = ".btl"

    def build(self) -> str | None:
        p = self.project
        if not p.back_analysis or not p.back_analysis.enabled:
            return None

        ba = p.back_analysis
        lines = []

        lines.append(f" {len(ba.parameters)}")
        for param in ba.parameters:
            lines.append(
                f" {param.material}  {param.property}  "
                f"{self._fmtg(param.initial)}  "
                f"{self._fmtg(param.bounds[0])}  {self._fmtg(param.bounds[1])}"
            )

        return "\n".join(lines) + "\n"
