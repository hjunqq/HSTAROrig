"""Generate the .ele (element connectivity) file.

Format: each line is
    elem_id  [spaces]  node1 node2 ... nodeN  group_id

Read by Global.f90:read_element() for each element in each group.
The elements must be ordered by group (group 1 elements first, then group 2, etc.).
"""

from __future__ import annotations

from .base import BaseGenerator


class EleGenerator(BaseGenerator):
    extension = ".ele"

    def build(self) -> str:
        if not self.mesh:
            return None

        lines = []
        # Elements must be in group order as read by global_data
        # First sort elements by group, preserving order within groups
        for eg in self.project.element_groups:
            gid = eg.group
            if gid in self.mesh.groups:
                for eid in self.mesh.groups[gid].element_ids:
                    elem = self.mesh.elements[eid]
                    nodes_str = "        ".join(f"{n}" for n in elem.nodes)
                    lines.append(
                        f"{eid:10d}                 {nodes_str}         {gid}"
                    )

        return "\n".join(lines) + "\n"
