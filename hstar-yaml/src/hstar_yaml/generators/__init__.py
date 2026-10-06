"""Generators package — generate all HSTAR input files from a project config."""

from __future__ import annotations

from pathlib import Path
import shutil

from ..schema import HstarProject
from ..mesh_reader import MeshData

from .inp import InpGenerator
from .glb import GlbGenerator
from .mat import MatGenerator
from .pre import PreGenerator
from .loa import LoaGenerator
from .man import ManGenerator
from .ele import EleGenerator
from .sol import SolGenerator
from .ftr import FtrGenerator
from .ini import IniGenerator
from .btl import BtlGenerator
from .obs import ObsGenerator


def generate_all(
    project: HstarProject,
    mesh: MeshData | None,
    output_dir: Path,
) -> list[Path]:
    """Generate all HSTAR input files and return list of generated paths."""
    generated: list[Path] = []

    generators = [
        InpGenerator(project, mesh),
        GlbGenerator(project, mesh),
        MatGenerator(project, mesh),
        PreGenerator(project, mesh),
        LoaGenerator(project, mesh),
        ManGenerator(project, mesh),
        EleGenerator(project, mesh),
        SolGenerator(project, mesh),
        FtrGenerator(project, mesh),
        IniGenerator(project, mesh),
        BtlGenerator(project, mesh),
        ObsGenerator(project, mesh),
    ]

    for gen in generators:
        try:
            path = gen.generate(output_dir)
            if path is not None:
                generated.append(path)
        except Exception as e:
            raise RuntimeError(
                f"Error generating {gen.__class__.__name__}: {e}"
            ) from e

    # Copy .msh file to output directory if not already there
    if project.mesh and project.mesh.file:
        msh_src = Path(project.mesh.file)
        if not msh_src.is_absolute():
            # Assume relative to output_dir or current dir
            pass
        msh_dst = output_dir / f"{project.problem_name}.msh"
        if msh_src.exists() and msh_src.resolve() != msh_dst.resolve():
            shutil.copy2(msh_src, msh_dst)
            generated.append(msh_dst)

    # Generate .cor file (coordinate file, extracted from mesh)
    cor_path = _generate_cor(project, mesh, output_dir)
    if cor_path:
        generated.append(cor_path)

    # Generate stub files that HSTAR expects to exist
    for ext in [".tem", ".ifs", ".nrt", ".aqu"]:
        stub_path = output_dir / f"{project.problem_name}{ext}"
        if not stub_path.exists():
            stub_path.write_text("", encoding="utf-8")
            generated.append(stub_path)

    return generated


def _generate_cor(
    project: HstarProject,
    mesh: MeshData | None,
    output_dir: Path,
) -> Path | None:
    """Generate .cor (coordinate) file from mesh data."""
    if not mesh or not mesh.nodes:
        return None

    path = output_dir / f"{project.problem_name}.cor"
    ndimn = project.mesh.dimension

    lines = []
    for nid in sorted(mesh.nodes.keys()):
        node = mesh.nodes[nid]
        if ndimn == 2:
            lines.append(f"{nid:8d}  {node.x:<16}  {node.y:<16}")
        else:
            lines.append(f"{nid:8d}  {node.x:<16}  {node.y:<16}  {node.z:<16}")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
