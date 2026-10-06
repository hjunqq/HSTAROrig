"""Generate the .man (management/analysis control) file.

Format (from Fem.f90 and Global.f90):
    For each increment block:
        Line 1: text " nincs,cdtest,earthquake_curve(1:ndimn)"
        Line 2: nincs cdtest eq_curve_x eq_curve_y [eq_curve_z]
        Line 3: nsteps dt istif icurve maxiter iprin outinc 0 0 0
        Line 4: tolerance values (11 values, can use shorthand like "11*1.e-5")
"""

from __future__ import annotations

from .base import BaseGenerator
from ..schema import AnalysisType


class ManGenerator(BaseGenerator):
    extension = ".man"

    def build(self) -> str:
        p = self.project
        a = p.analysis
        ndimn = p.mesh.dimension

        # Earthquake curve references (0 = none)
        eq_curves = [0] * ndimn
        if a.type == AnalysisType.DYNAMIC:
            for tc in p.time_curves:
                if tc.type.value == "SEISMIC":
                    dir_map = {"X": 0, "Y": 1, "Z": 2}
                    idx = dir_map.get(tc.direction.upper(), 0)
                    if idx < ndimn:
                        eq_curves[idx] = tc.id

        lines = []

        # If we have explicit increment blocks, use them
        if a.increments.blocks:
            for block in a.increments.blocks:
                lines.append(" nincs,cdtest,earthquake_curve(1:ndimn)")
                eq_str = "  ".join(str(c) for c in eq_curves)
                lines.append(f"  1  0  {eq_str}")

                nsteps = block.steps
                dt = block.dt
                istif = block.stiffness_update
                icurve = block.load_curve
                maxiter = block.max_iterations
                iprin = block.print_interval
                outinc = block.output_interval

                lines.append(
                    f"    {nsteps}      {dt}   {istif}   {icurve} {maxiter}   {iprin}   {outinc}   0 0 0 "
                )

                tol = block.tolerance.force
                lines.append(f" 11*{tol}")
        else:
            # Single block from legacy config
            lines.append(" nincs,cdtest,earthquake_curve(1:ndimn)")
            eq_str = "  ".join(str(c) for c in eq_curves)
            lines.append(f"  1  0  {eq_str}")

            if a.type == AnalysisType.DYNAMIC and a.dynamic:
                nsteps = int(a.dynamic.total_time / a.dynamic.dt) if a.dynamic.dt > 0 else a.increments.total_steps
                dt = a.dynamic.dt
            else:
                nsteps = a.increments.total_steps
                dt = a.increments.dt

            maxiter = a.increments.max_iterations
            lines.append(
                f"    {nsteps}     {dt}   1   1 {maxiter}  1   1   0 0 0 "
            )

            tol = a.increments.tolerance.force
            lines.append(f" 11*{tol}")

        return "\n".join(lines) + "\n"
