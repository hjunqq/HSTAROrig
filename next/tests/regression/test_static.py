"""Test static analysis cases"""
import subprocess, os, time

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'

test_dirs = [
    r'G:\BB\lyr\fix_boundary\changfang1',
    r'G:\BB\lyr\fix_boundary\chuanzha1',
    r'G:\BB\lyr\fix_boundary\youfei1',
    r'G:\BB\lyc\coarse_mesh\3cichuanzha\340initial',
]

for cwd in test_dirs:
    print(f'\n=== {os.path.basename(os.path.dirname(cwd))}/{os.path.basename(cwd)} ===')
    start = time.time()
    result = subprocess.run(
        [exe], cwd=cwd, capture_output=True, text=True, timeout=120, env=env,
    )
    elapsed = time.time() - start
    print(f'rc={result.returncode}, time={elapsed:.1f}s')
    if result.returncode != 0:
        stderr_lines = result.stderr.strip().splitlines()
        for line in stderr_lines[-3:]:
            print(f'  {line}')
    else:
        print('SUCCESS!')
        # Check chk file
        chk = os.path.join(cwd, '1.chk')
        if os.path.exists(chk):
            with open(chk, 'r') as f:
                lines = f.readlines()
            print(f'  1.chk: {len(lines)} lines')
