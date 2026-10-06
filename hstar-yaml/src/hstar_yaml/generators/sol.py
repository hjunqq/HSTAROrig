"""Generate the .sol (solver) file.

Format (from Solver.f90):
    For each block (nblks times, typically repeated):
        Line 1: text "mtype,ncpu,msglvl,"
        Line 2: mtype  ncpu  msglvl
        Line 3: text "isdefault  !0=use default;1=not use default"
        Line 4: isdefault (0 or 1)
        If isdefault == 1:
            Line 5: text "reducing_order,PreCGS,permutation,maxiter,out_of_core,eps_pivot,iparm11,iparm13"
            Line 6: values
"""

from __future__ import annotations

from .base import BaseGenerator
from ..schema import SolverType


class SolGenerator(BaseGenerator):
    extension = ".sol"

    def build(self) -> str:
        p = self.project
        s = p.solver
        nblks = max(len(p.construction_stages), 1)

        mtype = s.matrix_type
        ncpu = s.num_threads
        msglvl = s.message_level

        lines = []
        # Repeat for each block
        for _ in range(nblks):
            lines.append("mtype,ncpu,msglvl,")
            lines.append(f"  {mtype}    {ncpu}     {msglvl}")
            lines.append("isdefault  !0=use default;1=not use default")
            lines.append(f"    {s.isdefault}")

            if s.isdefault == 1:
                lines.append("reducing_order,PreCGS,permutation,maxiter,out_of_core,eps_pivot,iparm11,iparm13")
                lines.append(
                    f"      {s.reducing_order}          {s.PreCGS}        {s.permutation}         "
                    f"{s.maxiter}        {s.out_of_core}           {s.eps_pivot}        {s.iparm11}        {s.iparm13}"
                )

        return "\n".join(lines) + "\n"
