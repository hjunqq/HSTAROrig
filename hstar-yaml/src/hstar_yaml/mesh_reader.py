"""GiD .msh file reader for HSTAR.

GiD mesh format (as used by HSTAR):
    Line 1: header "NODES INFORMATION" (or similar)
    Node lines: node_id  x  y  [z]
    After nodes: "ELEMENT CONNECTIVITY" or similar header
    Element lines: elem_id  node1 node2 ... nodeN  group_id
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class MeshNode:
    id: int
    x: float
    y: float
    z: float = 0.0


@dataclass
class MeshElement:
    id: int
    nodes: list[int]
    group: int


@dataclass
class MeshGroup:
    """Tracks which elements belong to each group."""
    id: int
    element_ids: list[int] = field(default_factory=list)
    nelem: int = 0


@dataclass
class MeshData:
    """Parsed mesh data."""
    npoin: int = 0
    nelem: int = 0
    ndimn: int = 2
    nodes: dict[int, MeshNode] = field(default_factory=dict)
    elements: dict[int, MeshElement] = field(default_factory=dict)
    groups: dict[int, MeshGroup] = field(default_factory=dict)

    @property
    def max_node_id(self) -> int:
        return max(self.nodes.keys()) if self.nodes else 0

    @property
    def max_elem_id(self) -> int:
        return max(self.elements.keys()) if self.elements else 0


def read_gid_mesh(filepath: Path, ndimn: int = 2) -> MeshData:
    """Read a GiD .msh file and return parsed mesh data.

    The format is:
        NODES INFORMATION
        node_id  x  y  [z]
        ...
        (blank line or element section header)
        elem_id  [spaces]  n1 n2 n3 ... nN  group_id
        ...
    """
    mesh = MeshData(ndimn=ndimn)
    lines = Path(filepath).read_text(encoding="utf-8", errors="replace").splitlines()

    section = "unknown"
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue

        # Detect section headers
        upper = stripped.upper()
        if "NODE" in upper and ("INFO" in upper or "COORD" in upper):
            section = "nodes"
            continue
        if "ELEMENT" in upper and ("CONNECT" in upper or "INFO" in upper):
            section = "elements"
            continue
        if upper.startswith("END"):
            section = "unknown"
            continue

        # Parse based on section
        if section == "unknown":
            # Try to auto-detect: if first line looks like node data
            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    int(parts[0])
                    float(parts[1])
                    float(parts[2])
                    section = "nodes"
                except (ValueError, IndexError):
                    continue

        if section == "nodes":
            parts = stripped.split()
            if len(parts) < 3:
                continue
            try:
                node_id = int(parts[0])
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3]) if len(parts) > 3 else 0.0
                mesh.nodes[node_id] = MeshNode(id=node_id, x=x, y=y, z=z)
            except (ValueError, IndexError):
                # If we can't parse as a node, maybe we've moved to elements
                # Try to detect element section
                if len(parts) >= 4:
                    try:
                        all_ints = all(
                            parts[i].lstrip("-").isdigit()
                            for i in range(len(parts))
                        )
                        if all_ints and len(parts) >= 4:
                            section = "elements"
                            # Fall through to element parsing below
                        else:
                            continue
                    except (ValueError, IndexError):
                        continue
                else:
                    continue

        if section == "elements":
            parts = stripped.split()
            if len(parts) < 4:
                continue
            try:
                elem_id = int(parts[0])
                # Last value is group_id, middle values are node connectivity
                group_id = int(parts[-1])
                node_ids = [int(p) for p in parts[1:-1]]
                mesh.elements[elem_id] = MeshElement(
                    id=elem_id, nodes=node_ids, group=group_id
                )
                # Track group membership
                if group_id not in mesh.groups:
                    mesh.groups[group_id] = MeshGroup(id=group_id)
                mesh.groups[group_id].element_ids.append(elem_id)
            except (ValueError, IndexError):
                continue

    mesh.npoin = len(mesh.nodes)
    mesh.nelem = len(mesh.elements)

    # Update group element counts
    for g in mesh.groups.values():
        g.nelem = len(g.element_ids)

    return mesh


def read_cor_ele(cor_path: Path, ele_path: Path, ndimn: int = 3) -> MeshData:
    """Read HSTAR .cor + .ele files and return parsed mesh data.

    .cor format: node_id  x  y  [z]
    .ele format: elem_id  node1 node2 ... nodeN  group_id
    """
    mesh = MeshData(ndimn=ndimn)

    # Read coordinates
    cor_lines = Path(cor_path).read_text(encoding="utf-8", errors="replace").splitlines()
    for line in cor_lines:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 3:
            continue
        try:
            node_id = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3]) if len(parts) > 3 else 0.0
            mesh.nodes[node_id] = MeshNode(id=node_id, x=x, y=y, z=z)
        except (ValueError, IndexError):
            continue

    # Read elements
    ele_lines = Path(ele_path).read_text(encoding="utf-8", errors="replace").splitlines()
    for line in ele_lines:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            elem_id = int(parts[0])
            group_id = int(parts[-1])
            node_ids = [int(p) for p in parts[1:-1]]
            mesh.elements[elem_id] = MeshElement(
                id=elem_id, nodes=node_ids, group=group_id
            )
            if group_id not in mesh.groups:
                mesh.groups[group_id] = MeshGroup(id=group_id)
            mesh.groups[group_id].element_ids.append(elem_id)
        except (ValueError, IndexError):
            continue

    mesh.npoin = len(mesh.nodes)
    mesh.nelem = len(mesh.elements)

    for g in mesh.groups.values():
        g.nelem = len(g.element_ids)

    return mesh
