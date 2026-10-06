"""Test: does the existing hstar.exe in the source dir work?"""
import subprocess
import os

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

cwd = r'G:\BB\cube_contact\340(dizhen)'

# Test with original exe in the directory
exe = os.path.join(cwd, 'hstar.exe')
print(f'Testing: {exe}')
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
if result.returncode != 0:
    print('--- stderr (last 5 lines) ---')
    for line in result.stderr.strip().splitlines()[-5:]:
        print(f'  {line}')
else:
    print('SUCCESS')
