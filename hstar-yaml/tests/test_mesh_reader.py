"""Tests for GiD mesh reader."""

import pytest
import tempfile
from pathlib import Path

from hstar_yaml.mesh_reader import read_gid_mesh


SAMPLE_MSH_2D = """\
NODES INFORMATION
       1  0.0              0.0              0.0
       2  1.0              0.0              0.0
       3  1.0              1.0              0.0
       4  0.0              1.0              0.0
       5  0.5              0.0              0.0
       6  1.0              0.5              0.0
       7  0.5              1.0              0.0
       8  0.0              0.5              0.0
ELEMENT CONNECTIVITY
         1                  1         5         6         2         1
         2                  1         8         7         4         1
         3                  5         2         3         6         2
         4                  8         4         3         7         2
"""


class TestMeshReader:
    def test_read_basic_mesh(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write(SAMPLE_MSH_2D)
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert mesh.npoin == 8
        assert mesh.nelem == 4
        assert mesh.ndimn == 2

    def test_node_coordinates(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write(SAMPLE_MSH_2D)
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert mesh.nodes[1].x == 0.0
        assert mesh.nodes[1].y == 0.0
        assert mesh.nodes[3].x == 1.0
        assert mesh.nodes[3].y == 1.0

    def test_element_connectivity(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write(SAMPLE_MSH_2D)
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert mesh.elements[1].nodes == [1, 5, 6, 2]
        assert mesh.elements[1].group == 1

    def test_groups(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write(SAMPLE_MSH_2D)
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert 1 in mesh.groups
        assert 2 in mesh.groups
        assert mesh.groups[1].nelem == 2
        assert mesh.groups[2].nelem == 2

    def test_max_ids(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write(SAMPLE_MSH_2D)
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert mesh.max_node_id == 8
        assert mesh.max_elem_id == 4

    def test_empty_mesh(self):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".msh", delete=False, encoding="utf-8") as f:
            f.write("")
            f.flush()
            mesh = read_gid_mesh(Path(f.name), ndimn=2)

        assert mesh.npoin == 0
        assert mesh.nelem == 0
