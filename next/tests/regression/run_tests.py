"""
Regression test framework for HSTAR Phase 1 decomposition.

Verifies that the decomposed version produces identical output to the original.
Approach:
  1. Run reference (original) executable from the source directory
  2. Save output files
  3. Run test (decomposed) executable from the same source directory
  4. Compare output files with smart filtering for:
     - Timestamps (expected to differ)
     - Uninitialized memory values (denormalized near-zero floats)

Note: HSTAR must run from the directory containing its input files because
it reads input files using relative paths. Temp directory approach causes
access violations.

Usage:
  python run_tests.py [--ref-exe PATH] [--test-exe PATH] [--case CASE_NAME]
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent  # next/tests/regression -> project root

DEFAULT_REF_EXE = PROJECT_ROOT / 'next' / 'build_ref' / 'hstar_ref.exe'
DEFAULT_TEST_EXE = PROJECT_ROOT / 'next' / 'build' / 'hstar_phase1.exe'

# Test case definitions
TEST_CASES = {
    'cube_contact': {
        'source_dir': Path('G:/BB/cube_contact/340(dizhen)'),
        'problem_name': '1',
        'description': 'Contact analysis with earthquake loading - 40 nodes, 11 elements',
    },
}

# Output files to compare
OUTPUT_EXTENSIONS = ['.chk', '.flavia.res']

# ──────────────────────────────────────────────────────────────────────
# Smart comparison
# ──────────────────────────────────────────────────────────────────────

# Pattern for Fortran denormalized near-zero values like 0.36180-313 or 0.65977E-312
# These are uninitialized memory values, not meaningful computation results.
RE_DENORM = re.compile(r'[-+]?\d\.\d+[EeDd]?-(?:3[0-9]{2}|[4-9]\d{2})')

# Pattern for timestamp lines (date/time stamps vary between runs)
RE_TIMESTAMP = re.compile(
    r'(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)'
    r'|(?:\d{4}[/-]\d{2}[/-]\d{2})'
    r'|(?:\d{2}:\d{2}:\d{2})'
    r'|(?:Date\s*[&:=])'
    r'|(?:Time\s*[&:=])'
    r'|(?:Run\s+at)'
    r'|(?:elapsed)',
    re.IGNORECASE
)

# Pattern for Fortran-format numbers (including no-E notation like 0.12345-003)
RE_FORTRAN_NUM = re.compile(
    r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+|-\d{2,3}|\+\d{2,3})?'
)


def is_timestamp_line(line):
    """Check if a line contains timestamp information that's expected to differ."""
    return bool(RE_TIMESTAMP.search(line))


def normalize_denorms(line):
    """Replace denormalized near-zero values with a canonical '0.00000E+000'."""
    return RE_DENORM.sub('0.00000E+000', line)


def lines_equivalent(ref_line, test_line):
    """Check if two lines are equivalent after filtering known differences."""
    # Exact match
    if ref_line == test_line:
        return True

    # Both are timestamp lines - consider equivalent
    if is_timestamp_line(ref_line) and is_timestamp_line(test_line):
        return True

    # Normalize denormalized values and compare
    ref_norm = normalize_denorms(ref_line)
    test_norm = normalize_denorms(test_line)
    if ref_norm == test_norm:
        return True

    return False


def compare_output_file(ref_path, test_path, label=''):
    """Compare two output files with smart filtering.

    Returns (passed, summary_string).
    """
    if not ref_path.exists():
        return False, f'{label}: Reference file missing: {ref_path}'
    if not test_path.exists():
        return False, f'{label}: Test file missing: {test_path}'

    with open(ref_path, 'r', errors='replace') as f:
        ref_lines = f.readlines()
    with open(test_path, 'r', errors='replace') as f:
        test_lines = f.readlines()

    # Quick check
    if ref_lines == test_lines:
        return True, f'{label}: IDENTICAL ({len(ref_lines)} lines)'

    # Line count check
    if len(ref_lines) != len(test_lines):
        return False, (f'{label}: Line count differs '
                       f'(ref={len(ref_lines)}, test={len(test_lines)})')

    # Line-by-line smart comparison
    diffs = []
    n_timestamp_diffs = 0
    n_denorm_diffs = 0

    for i, (rl, tl) in enumerate(zip(ref_lines, test_lines), 1):
        if rl == tl:
            continue

        if is_timestamp_line(rl) or is_timestamp_line(tl):
            n_timestamp_diffs += 1
            continue

        ref_norm = normalize_denorms(rl)
        test_norm = normalize_denorms(tl)
        if ref_norm == test_norm:
            n_denorm_diffs += 1
            continue

        # Real difference
        diffs.append((i, rl.rstrip(), tl.rstrip()))

    if not diffs:
        parts = [f'{label}: EQUIVALENT ({len(ref_lines)} lines)']
        if n_timestamp_diffs:
            parts.append(f'{n_timestamp_diffs} timestamp diffs (filtered)')
        if n_denorm_diffs:
            parts.append(f'{n_denorm_diffs} denorm diffs (filtered)')
        return True, ', '.join(parts)
    else:
        parts = [f'{label}: {len(diffs)} REAL DIFFERENCES found']
        if n_timestamp_diffs:
            parts.append(f'({n_timestamp_diffs} timestamp diffs filtered)')
        if n_denorm_diffs:
            parts.append(f'({n_denorm_diffs} denorm diffs filtered)')
        summary = ', '.join(parts)
        detail_lines = []
        for lineno, ref_l, test_l in diffs[:15]:
            detail_lines.append(f'  Line {lineno}:')
            detail_lines.append(f'    ref:  {ref_l}')
            detail_lines.append(f'    test: {test_l}')
        return False, summary + '\n' + '\n'.join(detail_lines)


# ──────────────────────────────────────────────────────────────────────
# Execution helpers
# ──────────────────────────────────────────────────────────────────────

def get_intel_env():
    """Create environment with Intel runtime DLL paths."""
    env = os.environ.copy()
    intel_paths = [
        r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
        r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
        r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
    ]
    existing_path = env.get('PATH', '')
    env['PATH'] = ';'.join(intel_paths) + ';' + existing_path
    return env


def run_exe(exe_path, work_dir, env, timeout=300):
    """Run the HSTAR executable. Returns (returncode, stdout, stderr, elapsed)."""
    start = time.time()
    try:
        result = subprocess.run(
            [str(exe_path)],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
        elapsed = time.time() - start
        return result.returncode, result.stdout, result.stderr, elapsed
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start
        return -1, '', f'TIMEOUT after {timeout}s', elapsed


def collect_outputs(source_dir, problem_name):
    """Collect output file paths that exist."""
    outputs = {}
    for ext in OUTPUT_EXTENSIONS:
        p = source_dir / (problem_name + ext)
        if p.exists():
            outputs[ext] = p
    return outputs


def remove_outputs(source_dir, problem_name):
    """Remove output files from the source directory."""
    for ext in OUTPUT_EXTENSIONS:
        p = source_dir / (problem_name + ext)
        if p.exists():
            p.unlink()


def save_outputs(source_dir, problem_name, dest_dir):
    """Copy output files to a destination directory."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    saved = {}
    for ext in OUTPUT_EXTENSIONS:
        p = source_dir / (problem_name + ext)
        if p.exists():
            dest = dest_dir / (problem_name + ext)
            shutil.copy2(p, dest)
            saved[ext] = dest
    return saved


# ──────────────────────────────────────────────────────────────────────
# Main test runner
# ──────────────────────────────────────────────────────────────────────

def run_test_case(case_name, case_def, ref_exe, test_exe):
    """Run a single test case. Returns (passed, details)."""
    print(f'\n  Test: {case_name}')
    print(f'    {case_def["description"]}')

    source_dir = case_def['source_dir']
    problem_name = case_def['problem_name']

    if not source_dir.exists():
        return False, f'Source directory not found: {source_dir}'

    # Create temp directory for saving outputs
    results_dir = Path(tempfile.mkdtemp(prefix=f'hstar_results_{case_name}_'))
    ref_results = results_dir / 'ref'
    test_results = results_dir / 'test'

    env = get_intel_env()
    all_pass = False  # Will be set to True only if all comparisons pass

    # Backup existing output files (if any) so we can restore them later
    backup_dir = results_dir / 'backup'
    backup_dir.mkdir()
    existing_outputs = collect_outputs(source_dir, problem_name)
    for ext, p in existing_outputs.items():
        shutil.copy2(p, backup_dir / p.name)

    try:
        # ── Step 1: Run reference executable ──
        remove_outputs(source_dir, problem_name)
        print(f'    Running reference: {ref_exe.name} ...')
        rc_ref, out_ref, err_ref, t_ref = run_exe(ref_exe, source_dir, env)
        print(f'    Reference: rc={rc_ref}, time={t_ref:.1f}s')

        if rc_ref != 0:
            (results_dir / 'ref_stdout.txt').write_text(out_ref)
            (results_dir / 'ref_stderr.txt').write_text(err_ref)
            return False, f'Reference exe failed (rc={rc_ref}). Logs: {results_dir}'

        ref_saved = save_outputs(source_dir, problem_name, ref_results)
        print(f'    Reference outputs saved: {list(ref_saved.keys())}')

        # ── Step 2: Run test executable ──
        remove_outputs(source_dir, problem_name)
        print(f'    Running test: {test_exe.name} ...')
        rc_test, out_test, err_test, t_test = run_exe(test_exe, source_dir, env)
        print(f'    Test: rc={rc_test}, time={t_test:.1f}s')

        if rc_test != 0:
            (results_dir / 'test_stdout.txt').write_text(out_test)
            (results_dir / 'test_stderr.txt').write_text(err_test)
            return False, f'Test exe failed (rc={rc_test}). Logs: {results_dir}'

        test_saved = save_outputs(source_dir, problem_name, test_results)
        print(f'    Test outputs saved: {list(test_saved.keys())}')

        # ── Step 3: Compare outputs ──
        results = []
        all_pass = True

        for ext in OUTPUT_EXTENSIONS:
            ref_f = ref_results / (problem_name + ext)
            test_f = test_results / (problem_name + ext)

            if not ref_f.exists() and not test_f.exists():
                results.append(f'{ext}: not produced by either run (skipped)')
                continue

            match, detail = compare_output_file(ref_f, test_f, label=ext)
            results.append(detail)
            if not match:
                all_pass = False

        detail_str = '\n    '.join(results)
        return all_pass, detail_str

    finally:
        # Restore original output files
        remove_outputs(source_dir, problem_name)
        for p in backup_dir.iterdir():
            shutil.copy2(p, source_dir / p.name)

        # Clean up results dir if all passed
        if all_pass:
            shutil.rmtree(results_dir, ignore_errors=True)
            print(f'    Cleaned up results dir')
        else:
            print(f'    Results preserved at: {results_dir}')


def main():
    parser = argparse.ArgumentParser(description='HSTAR Phase 1 Regression Tests')
    parser.add_argument('--ref-exe', type=Path, default=DEFAULT_REF_EXE,
                        help='Path to reference (original) executable')
    parser.add_argument('--test-exe', type=Path, default=DEFAULT_TEST_EXE,
                        help='Path to test (decomposed) executable')
    parser.add_argument('--case', type=str, default=None,
                        help='Run specific test case (default: all)')
    args = parser.parse_args()

    # Validate executables
    if not args.ref_exe.exists():
        print(f'ERROR: Reference executable not found: {args.ref_exe}')
        sys.exit(1)
    if not args.test_exe.exists():
        print(f'ERROR: Test executable not found: {args.test_exe}')
        sys.exit(1)

    print('=' * 60)
    print('HSTAR Phase 1 Regression Tests')
    print('=' * 60)
    print(f'Reference exe: {args.ref_exe}')
    print(f'Test exe:      {args.test_exe}')

    cases = TEST_CASES
    if args.case:
        if args.case not in cases:
            print(f'ERROR: Unknown case "{args.case}". Available: {list(cases.keys())}')
            sys.exit(1)
        cases = {args.case: cases[args.case]}

    passed = 0
    failed = 0
    summaries = []

    for name, case_def in cases.items():
        ok, detail = run_test_case(name, case_def, args.ref_exe, args.test_exe)
        if ok:
            passed += 1
            summaries.append(f'  PASS: {name}')
        else:
            failed += 1
            summaries.append(f'  FAIL: {name}')
        print(f'    Result: {"PASS" if ok else "FAIL"}')
        print(f'    {detail}')

    print('\n' + '=' * 60)
    print('Summary')
    print('=' * 60)
    for s in summaries:
        print(s)
    print(f'\n  Total: {passed + failed}, Passed: {passed}, Failed: {failed}')
    print('=' * 60)

    sys.exit(0 if failed == 0 else 1)


if __name__ == '__main__':
    main()
