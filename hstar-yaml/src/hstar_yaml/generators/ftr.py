"""Generate the .ftr (surface force/safety factor) file.

Format (from Global.f90 L834+):
    Line 1: text
    Line 2: nforce ngaps nforce_gaps nsafety_gaps
    If nforce != 0:
        Line: text "nforce_appear"
        Line: nforce_appear values
        Line: text
        For each force surface:
            Line: text
            Line: lgroup neface node_face nliste
            Line: text
            Line: group list
            For each face element:
                Line: elem_id  face_node1 face_node2 [face_node3...]
"""

from __future__ import annotations

from .base import BaseGenerator


class FtrGenerator(BaseGenerator):
    extension = ".ftr"

    def build(self) -> str:
        p = self.project
        sf = p.surface_forces

        lines = []
        lines.append(" nforce,ngaps,nforce_gaps,nsafety_gaps")

        if sf and sf.forces:
            nforce = len(sf.forces)
            lines.append(f"       {nforce}       {sf.safety_gaps}       0       0")

            lines.append("nforce_appear !1-- for saftyfactor 2-- for internal force 3-- for both")
            appear_vals = " ".join(str(f.appear_type) for f in sf.forces)
            lines.append(appear_vals)

            lines.append(" 1:nforce")
            for i, force in enumerate(sf.forces, 1):
                lines.append(f" lgroup,neface,node_face,nliste !iforce=           {i}")
                neface = len(force.elements)
                lines.append(
                    f"  {force.group}       {neface}     {force.face_nodes}           1 "
                )
                lines.append(" surface_force(iforce)%list")
                lines.append(f"           {force.group}")

                for eid in force.elements:
                    lines.append(f"      {eid}")
        else:
            lines.append("       0       0       0       0")

        return "\n".join(lines) + "\n"
