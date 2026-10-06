"""Generate the .obs (observation) file for back-analysis.

Format:
    Line 1: nobs
    For each observation:
        Line: node  dof  weight
    Observation data file reference
"""

from __future__ import annotations

from .base import BaseGenerator
from ..schema import DofName


DOF_INDEX = {
    DofName.UX: 1,
    DofName.UY: 2,
    DofName.UZ: 3,
    DofName.PW: 8,
    DofName.T: 10,
}


class ObsGenerator(BaseGenerator):
    extension = ".obs"

    def build(self) -> str | None:
        p = self.project
        if not p.back_analysis or not p.back_analysis.enabled:
            return None
        if not p.back_analysis.observations:
            return None

        obs = p.back_analysis.observations
        points = obs.get("points", [])

        lines = []
        lines.append(f" {len(points)}")
        for pt in points:
            node = pt.get("node", 0)
            dof = pt.get("dof", "UY")
            weight = pt.get("weight", 1.0)
            dof_idx = DOF_INDEX.get(DofName(dof), 2)
            lines.append(f" {node}  {dof_idx}  {self._fmtg(weight)}")

        data_file = obs.get("file", "")
        if data_file:
            lines.append(f" {data_file}")

        return "\n".join(lines) + "\n"
