"""End-to-end comparison: generate from 340initial.yaml and compare with original files."""

import sys
import os
import tempfile
from pathlib import Path

# Fix Windows console encoding
if sys.platform == "win32":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from hstar_yaml.schema import HstarProject
from hstar_yaml.mesh_reader import read_cor_ele
from hstar_yaml.generators.glb import GlbGenerator
from hstar_yaml.generators.mat import MatGenerator
from hstar_yaml.generators.sol import SolGenerator
from hstar_yaml.generators.man import ManGenerator

import yaml

CASE_DIR = Path(r"G:\BB\lyc\coarse_mesh\3cichuanzha\340initial")
YAML_FILE = Path(__file__).parent / "340initial.yaml"


def read_original(name: str) -> str:
    return (CASE_DIR / name).read_text(encoding="utf-8", errors="replace")


def compare_sections(name: str, original: str, generated: str):
    """Compare key sections of two files."""
    orig_lines = original.strip().splitlines()
    gen_lines = generated.strip().splitlines()

    diffs = []
    max_lines = max(len(orig_lines), len(gen_lines))
    for i in range(max_lines):
        ol = orig_lines[i].strip() if i < len(orig_lines) else "<MISSING>"
        gl = gen_lines[i].strip() if i < len(gen_lines) else "<MISSING>"
        if ol != gl:
            diffs.append((i + 1, ol, gl))

    if diffs:
        print(f"\n{'='*60}")
        print(f"DIFF: {name} — {len(diffs)} line(s) differ out of {max_lines}")
        print(f"{'='*60}")
        for lineno, ol, gl in diffs[:30]:
            print(f"  Line {lineno}:")
            print(f"    ORIG: {ol[:140]}")
            print(f"    GEN:  {gl[:140]}")
        if len(diffs) > 30:
            print(f"  ... and {len(diffs) - 30} more differences")
    else:
        print(f"  OK: {name} — matches perfectly ({max_lines} lines)")

    return len(diffs)


def main():
    print("340initial End-to-End Comparison Test")
    print("=" * 60)

    # Load YAML
    with open(YAML_FILE, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    project = HstarProject.model_validate(raw)

    # Read mesh from .cor + .ele
    cor_path = CASE_DIR / "1.cor"
    ele_path = CASE_DIR / "1.ele"
    mesh = read_cor_ele(cor_path, ele_path, ndimn=3)
    print(f"Mesh: {mesh.npoin} nodes, {mesh.nelem} elements, {len(mesh.groups)} groups")
    for gid, grp in sorted(mesh.groups.items()):
        print(f"  Group {gid}: {grp.nelem} elements")

    total_diffs = 0

    # Compare .glb
    print("\n--- Comparing 1.glb ---")
    gen = GlbGenerator(project, mesh)
    generated_glb = gen.build()
    total_diffs += compare_sections("1.glb", read_original("1.glb"), generated_glb)

    # Compare .mat
    print("\n--- Comparing 1.mat ---")
    gen = MatGenerator(project, mesh)
    generated_mat = gen.build()
    total_diffs += compare_sections("1.mat", read_original("1.mat"), generated_mat)

    # Compare .sol
    print("\n--- Comparing 1.sol ---")
    gen = SolGenerator(project, mesh)
    generated_sol = gen.build()
    total_diffs += compare_sections("1.sol", read_original("1.sol"), generated_sol)

    # Compare .man
    print("\n--- Comparing 1.man ---")
    gen = ManGenerator(project, mesh)
    generated_man = gen.build()
    total_diffs += compare_sections("1.man", read_original("1.man"), generated_man)

    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: {total_diffs} total line differences across all files")
    if total_diffs == 0:
        print("ALL FILES MATCH!")
    else:
        print("Some differences remain — review above.")


if __name__ == "__main__":
    main()
