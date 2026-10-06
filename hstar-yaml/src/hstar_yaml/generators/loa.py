"""Generate the .loa (load) file.

Format (from Load.f90:external_load_1()):
    Line 1: text (title)
    Line 2: ntcurve (number of time curves)
    For each time curve:
        Line: ncurve_id, type_str, 0, npoints
        For LINEAR/TABULAR:
            Lines: t_i  val_i pairs
        For SEISMIC:
            Line: nlines, 'SEISMIC', 0, ndata_lines
            Line: dt  start_val  total_time  scale_factor
            Lines: acceleration values (one per line)
    After curves:
        Point loads section
        Distributed loads section
"""

from __future__ import annotations

from pathlib import Path

from .base import BaseGenerator


class LoaGenerator(BaseGenerator):
    extension = ".loa"

    def build(self) -> str:
        p = self.project

        lines = []
        # Header
        lines.append(f"({p.problem_name}.loa)")

        # Number of time curves
        ntcurve = len(p.time_curves)
        lines.append(f"         {ntcurve}")

        # Time curves
        for tc in p.time_curves:
            if tc.type.value == "LINEAR" or tc.type.value == "TABULAR":
                npoints = len(tc.points)
                lines.append(f"    {tc.id}   LINEAR                  0              {npoints}")
                for pt in tc.points:
                    lines.append(f"   {pt[0]}   {pt[1]}")
            elif tc.type.value == "SEISMIC":
                # Seismic curve: read from file or inline
                nlines = tc.nlines
                lines.append(
                    f"{tc.id * 1000 + 1},'SEISMIC',0 ,{nlines} (nline)  !   x-a   ntime for curve {tc.id}"
                )
                lines.append(
                    f"{tc.dt}  0.0  {tc.total_time}  {tc.scale}  ! time"
                )
                # If file provided, include note
                if tc.file:
                    lines.append(f"! Seismic data from file: {tc.file}")
                    lines.append("! Include seismic acceleration data below, one value per line")

        # Gravity load (if any) - typically handled through body forces in HSTAR
        # not directly in the .loa file, but through NGRAV flag in .glb

        return "\n".join(lines) + "\n"
