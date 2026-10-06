"""Test smaller static cases with longer timeout"""
import subprocess, os, time

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'

# Try initial (static) cases - likely smaller
test_dirs = [
    r'G:\BB\lyc\coarse_mesh\3cichuanzha\340initial',
    r'G:\BB\lyc\sparse_mesh\340initial',
    r'G:\BB\lyc\contactjisuanwenjian\340initial',
]

for cwd in test_dirs:
    if not os.path.exists(cwd):
        print(f'SKIP (not found): {cwd}')
        continue
    # Check inp format
    inp = os.path.join(cwd, 'inp')
    with open(inp, 'r') as f:
        content = f.read()
    print(f'\n=== {os.path.basename(os.path.dirname(cwd))}/{os.path.basename(cwd)} ===')
    print(f'  inp: {content.strip()[:100]}')

    # Check msh file size
    msh = os.path.join(cwd, '1.msh')
    if os.path.exists(msh):
        print(f'  1.msh size: {os.path.getsize(msh)} bytes')

    start = time.time()
    result = subprocess.run(
        [exe], cwd=cwd, capture_output=True, text=True, timeout=300, env=env,
    )
    elapsed = time.time() - start
    print(f'  rc={result.returncode}, time={elapsed:.1f}s')
    if result.returncode != 0:
        stderr_lines = result.stderr.strip().splitlines()
        for line in stderr_lines[-5:]:
            print(f'    {line}')
    else:
        print('  SUCCESS!')
        chk = os.path.join(cwd, '1.chk')
        if os.path.exists(chk):
            with open(chk, 'r') as f:
                lines = f.readlines()
            print(f'  1.chk: {len(lines)} lines')
