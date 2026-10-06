"""Generate the .ini (initial stress) file.

Format:
    Line 1: text "STRESS" or "NONE"
    Line 2: ngroup_init  kinit_type
    Line 3: igroup  flag  extra_params...
    For each element gauss point:
        Line: elem_id  gauss_id  stress_components...
"""

from __future__ import annotations

from .base import BaseGenerator


class IniGenerator(BaseGenerator):
    extension = ".ini"

    def build(self) -> str:
        # Generate a minimal .ini file (no initial stress)
        # User would typically generate this from a previous analysis
        p = self.project
        ngroup = len(p.element_groups)
        ndimn = p.mesh.dimension

        lines = []
        lines.append("                              STRESS")
        lines.append(f"           {ngroup}           1")

        # For each group, write zero initial stress
        for eg in p.element_groups:
            lines.append(
                f"           {eg.group}           0  "
                f"0.000000000000000E+000  0.000000000000000E+000"
            )

        return "\n".join(lines) + "\n"
