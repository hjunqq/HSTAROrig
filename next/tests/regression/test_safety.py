"""Test chuanzha-safety static case"""
import subprocess, os, time

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

ref_exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'
test_exe = r'G:\BB\hstarYLOrig\next\build\hstar_phase1.exe'

cwd = r'G:\BB\lyr\fix_boundary\chuanzha-safety'

for label, exe in [('REF', ref_exe), ('PHASE1', test_exe)]:
    print(f'\n=== {label}: {os.path.basename(exe)} ===')
    start = time.time()
    result = subprocess.run(
        [exe], cwd=cwd, capture_output=True, text=True, timeout=600, env=env,
    )
    elapsed = time.time() - start
    print(f'  rc={result.returncode}, time={elapsed:.1f}s')
    if result.returncode != 0:
        for line in result.stderr.strip().splitlines()[-5:]:
            print(f'    {line}')
    else:
        print('  SUCCESS!')
        chk = os.path.join(cwd, '1.chk')
        if os.path.exists(chk):
            with open(chk, 'r') as f:
                lines = f.readlines()
            print(f'  1.chk: {len(lines)} lines, {os.path.getsize(chk)} bytes')
