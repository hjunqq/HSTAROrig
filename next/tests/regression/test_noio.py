"""Test running HSTAR without output capture"""
import subprocess
import os

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'
cwd = r'G:\BB\cube_contact\340(dizhen)'

# Test: no output capture, redirect to files
print('Running without capture_output...')
with open(os.path.join(cwd, '_test_stdout.txt'), 'w') as fout, \
     open(os.path.join(cwd, '_test_stderr.txt'), 'w') as ferr:
    result = subprocess.run(
        [exe], cwd=cwd, stdout=fout, stderr=ferr, timeout=120,
        env=env,
    )
print(f'rc={result.returncode}')

# Read results
with open(os.path.join(cwd, '_test_stderr.txt'), 'r') as f:
    stderr = f.read().strip()
if stderr:
    print('stderr:')
    for line in stderr.splitlines()[-5:]:
        print(f'  {line}')

# Check if chk file was updated
import time
chk_path = os.path.join(cwd, '1.chk')
if os.path.exists(chk_path):
    mtime = os.path.getmtime(chk_path)
    age = time.time() - mtime
    print(f'1.chk age: {age:.0f}s')
    with open(chk_path, 'r') as f:
        lines = f.readlines()
    print(f'1.chk lines: {len(lines)}')

# Clean up test files
os.remove(os.path.join(cwd, '_test_stdout.txt'))
os.remove(os.path.join(cwd, '_test_stderr.txt'))
