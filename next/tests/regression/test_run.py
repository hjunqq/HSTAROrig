"""Quick test: can hstar_ref.exe run via subprocess from source dir?"""
import subprocess
import os
import sys

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'
cwd = r'G:\BB\cube_contact\340(dizhen)'

print(f'Running {exe}')
print(f'  from {cwd}')

result = subprocess.run(
    [exe],
    cwd=cwd,
    capture_output=True,
    text=True,
    timeout=120,
    env=env,
)
print(f'rc={result.returncode}')
print(f'stdout lines: {len(result.stdout.splitlines())}')
print(f'stderr lines: {len(result.stderr.splitlines())}')
if result.returncode != 0:
    print('--- stderr (last 10 lines) ---')
    for line in result.stderr.strip().splitlines()[-10:]:
        print(f'  {line}')
    print('--- stdout (last 10 lines) ---')
    for line in result.stdout.strip().splitlines()[-10:]:
        print(f'  {line}')
else:
    print('SUCCESS')
    print('--- stdout (last 5 lines) ---')
    for line in result.stdout.strip().splitlines()[-5:]:
        print(f'  {line}')
