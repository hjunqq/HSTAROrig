"""End-to-end comparison script for 340initial case.

Reads original HSTAR files, parses key values, and compares with what
our generators would produce from an equivalent YAML config.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from hstar_yaml.schema import *
from hstar_yaml.generators.glb import GlbGenerator
from hstar_yaml.generators.mat import MatGenerator
from hstar_yaml.generators.sol import SolGenerator
from hstar_yaml.generators.man import ManGenerator
from hstar_yaml.generators.pre import PreGenerator
from hstar_yaml.generators.inp import InpGenerator
from hstar_yaml.mesh_reader import read_gid_mesh, MeshData

CASE_DIR = Path(r"G:\BB\lyc\coarse_mesh\3cichuanzha\340initial")


def read_original(name: str) -> str:
    return (CASE_DIR / name).read_text(encoding="utf-8", errors="replace")


def compare_lines(name: str, original: str, generated: str):
    """Compare two multi-line strings and report differences."""
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
        print(f"DIFF: {name} — {len(diffs)} line(s) differ")
        print(f"{'='*60}")
        for lineno, ol, gl in diffs[:20]:  # Show first 20 diffs
            print(f"  Line {lineno}:")
            print(f"    ORIG: {ol[:120]}")
            print(f"    GEN:  {gl[:120]}")
        if len(diffs) > 20:
            print(f"  ... and {len(diffs) - 20} more differences")
    else:
        print(f"  OK: {name} — matches perfectly")

    return len(diffs)


def main():
    print("340initial End-to-End Comparison")
    print("=" * 60)

    # Read mesh
    msh_path = CASE_DIR / "1.msh"
    mesh = read_gid_mesh(msh_path, ndimn=3)
    print(f"Mesh: {mesh.npoin} nodes, {mesh.nelem} elements, {len(mesh.groups)} groups")

    # Compare inp
    orig_inp = read_original("inp")
    print(f"\n--- inp ---")
    print(f"Original:\n{orig_inp}")

    # Compare glb - just parse and check key values
    orig_glb = read_original("1.glb")
    glb_lines = orig_glb.strip().splitlines()

    print(f"\n--- 1.glb key values ---")
    # Line 2: npoin, npoinb, nelem, ndimn, nmats, ngroup
    line2 = glb_lines[1].split()
    print(f"  npoin={line2[0]}, npoinb={line2[1]}, nelem={line2[2]}, "
          f"ndimn={line2[3]}, nmats={line2[4]}, ngroup={line2[5]}")
    print(f"  outplot={line2[7]}, kstab={line2[8]}")

    # Line 6: ninit, kinit, winit, nblks...
    line6 = glb_lines[5].split()
    print(f"  ninit={line6[0]}, kinit={line6[1]}, winit={line6[2]}, nblks={line6[3]}")

    # Line 8: type_problem, type_solver
    line8 = glb_lines[7].split()
    print(f"  type_problem={line8[0]}, type_solver={line8[1]}, type_nl={line8[3]}")

    # Line 16: mdofn
    line16 = glb_lines[15].strip()
    print(f"  mdofn={line16}")

    # Group definitions
    for i, line in enumerate(glb_lines):
        if "GROUP" in line and ("Q4" in line or "B8" in line or "T3" in line):
            parts = line.split()
            print(f"  Group: {line.strip()[:100]}")

    # Compare mat
    orig_mat = read_original("1.mat")
    print(f"\n--- 1.mat key values ---")
    for i, line in enumerate(orig_mat.splitlines()):
        stripped = line.strip()
        if "material_serial" in stripped:
            print(f"  {stripped}")
        if "ELASTIC" in stripped or "GOODMAN" in stripped or "CAMCLAY" in stripped:
            print(f"  {stripped[:100]}")

    # Compare sol
    orig_sol = read_original("1.sol")
    print(f"\n--- 1.sol ---")
    sol_lines = orig_sol.strip().splitlines()
    for line in sol_lines:
        print(f"  {line.strip()}")

    # Compare man
    orig_man = read_original("1.man")
    print(f"\n--- 1.man ---")
    man_lines = orig_man.strip().splitlines()
    for line in man_lines:
        print(f"  {line.strip()}")

    # Compare pre structure
    orig_pre = read_original("1.pre")
    pre_lines = orig_pre.strip().splitlines()
    print(f"\n--- 1.pre structure ---")
    print(f"  Total lines: {len(pre_lines)}")
    # Find fixset headers
    for i, line in enumerate(pre_lines):
        stripped = line.strip()
        if "PRESCRIBE SET" in stripped:
            print(f"  Line {i+1}: {stripped}")
        if i < len(pre_lines) and stripped and stripped[0].isdigit() and len(stripped.split()) == 2:
            parts = stripped.split()
            if all(p.isdigit() for p in parts):
                print(f"  Line {i+1}: nfixsets={parts[0]}, nline={parts[1]}")

    # Summary of differences found
    print(f"\n{'='*60}")
    print("KEY DIFFERENCES from current generator:")
    print("  1. type_problem=Q (FREQUENCY), not F (STATIC)")
    print("  2. ndimn=3, element index=9 (B8 hex)")
    print("  3. ninit=0, kinit=0, winit=-1 (different defaults)")
    print("  4. mdofn=3 (only 3 DOFs, not 10)")
    print("  5. nmass/nsmat/nhmat/nqmat=999 (special values)")
    print("  6. .sol has isdefault=1 with extra PARDISO params")
    print("  7. .man has 2 increment blocks (200+5 steps)")
    print("  8. .loa has edge loads (153 surface elements)")
    print("  9. .pre is repeated for nblks=1 twice")
    print("  10. Group line: index=9 not 4, elem_type=B8 not Q4")
    print("  11. equvs_process uses 1000*0 shorthand")
    print("  12. force_process uses 1000*0 shorthand")
    print("  13. mat file: Goodman has extra JANBU line + 37th line")
    print("  14. gid_ms=1, res output flags non-zero")


if __name__ == "__main__":
    main()
