"""Generate the .pre (prescribed boundary conditions) file.

Format (from Prescrib.f90:prescrib_set()):
    For each block (nblks times):
        Line: text header "PRESCRIBE SET--NFIXSETS, nblks"
        Line: nfixsets  nline
        For each fixset:
            Line: ifixvar nfixnods itcurve tfixvar outfix jfixvar gamawx nextr
                  (if type_ABC=='MIF': ifixvar ifixvar0 nfixnods ...)
            Line: node_list (nfixnods integers)
            Line: value_list (nfixnods reals, can use N*val shorthand)
"""

from __future__ import annotations

from .base import BaseGenerator
from ..schema import DofName, BcType


# Map DOF names to HSTAR ifixvar index (1-based)
DOF_TO_FIXVAR = {
    DofName.UX: 1,
    DofName.UY: 2,
    DofName.UZ: 3,
    DofName.RX: 4,
    DofName.RY: 5,
    DofName.RZ: 6,
    DofName.PW: 8,
    DofName.T: 10,
}


class PreGenerator(BaseGenerator):
    extension = ".pre"

    def build(self) -> str:
        p = self.project
        nblks = max(len(p.construction_stages), 1)

        # Group boundary conditions by DOF (ifixvar)
        # Each unique ifixvar becomes one fixset
        fixsets = self._build_fixsets()
        nfixsets = len(fixsets)

        # Compute total lines per block (for nline field)
        # Each fixset has: 1 param line + 1 node list line + 1 value line = 3 lines
        nline = nfixsets * 3

        lines = []

        # Repeat for each block
        for _iblk in range(nblks):
            lines.append(
                f"                                    PRESCRIBE SET--NFIXSETS,    {nblks}"
            )
            lines.append(f"       {nfixsets}      {nline}")

            for fs in fixsets:
                ifixvar = fs["ifixvar"]
                nodes = fs["nodes"]
                values = fs["values"]
                nfixnods = len(nodes)
                itcurve = fs.get("time_curve", 0)
                tfixvar = 0  # 0=displacement, 1=velocity, 2=acceleration
                outfix = 0
                jfixvar = 0
                gamawx = 0.0
                nextr = 0

                # ifixvar line
                lines.append(
                    f"         {ifixvar}        {nfixnods}         {itcurve}         {tfixvar}         {outfix}         {jfixvar}   {self._fmtg(gamawx)}         {nextr}"
                )

                # Node list
                node_str = " ".join(str(n) for n in nodes)
                lines.append(node_str)

                # Value list - use shorthand if all same
                if len(set(values)) == 1:
                    lines.append(f"{nfixnods}*{self._fmtg(values[0])}")
                else:
                    val_str = " ".join(self._fmtg(v) for v in values)
                    lines.append(val_str)

                lines.append("")  # blank line after each fixset

        return "\n".join(lines) + "\n"

    def _build_fixsets(self) -> list[dict]:
        """Build fixsets from boundary conditions, one per unique (dof, type) combo."""
        fixsets = []

        for bc in self.project.boundary_conditions:
            # Expand node_range if used
            nodes = list(bc.nodes)
            if bc.node_range:
                start = bc.node_range.get("start", 1)
                end = bc.node_range.get("end", 1)
                step = bc.node_range.get("step", 1)
                nodes.extend(range(start, end + 1, step))

            if not nodes:
                continue

            for dof in bc.dofs:
                ifixvar = DOF_TO_FIXVAR.get(dof, 1)
                values = [bc.value] * len(nodes)
                tc = bc.time_curve if bc.time_curve is not None else 0

                fixsets.append({
                    "ifixvar": ifixvar,
                    "nodes": nodes,
                    "values": values,
                    "time_curve": tc,
                })

        return fixsets
