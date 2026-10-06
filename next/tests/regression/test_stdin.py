"""Test different stdin approaches for running HSTAR"""
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

# Test 1: stdin=DEVNULL
print('=== Test 1: stdin=DEVNULL ===')
result = subprocess.run(
    [exe], cwd=cwd, capture_output=True, text=True, timeout=120,
    env=env, stdin=subprocess.DEVNULL,
)
print(f'rc={result.returncode}')
if result.returncode != 0:
    for line in result.stderr.strip().splitlines()[-3:]:
        print(f'  {line}')

# Test 2: stdin from inp file
print('\n=== Test 2: stdin from inp file ===')
inp_path = os.path.join(cwd, 'inp')
with open(inp_path, 'r') as inp_f:
    result = subprocess.run(
        [exe], cwd=cwd, capture_output=True, text=True, timeout=120,
        env=env, stdin=inp_f,
    )
print(f'rc={result.returncode}')
if result.returncode != 0:
    for line in result.stderr.strip().splitlines()[-3:]:
        print(f'  {line}')
else:
    print('SUCCESS!')
    for line in result.stdout.strip().splitlines()[-5:]:
        print(f'  {line}')

# Test 3: using cmd /c with redirect
print('\n=== Test 3: cmd /c ===')
result = subprocess.run(
    f'cmd /c "{exe}"',
    cwd=cwd, capture_output=True, text=True, timeout=120,
    env=env, shell=True,
)
print(f'rc={result.returncode}')
if result.returncode != 0:
    for line in result.stderr.strip().splitlines()[-3:]:
        print(f'  {line}')
else:
    print('SUCCESS!')
    for line in result.stdout.strip().splitlines()[-5:]:
        print(f'  {line}')
