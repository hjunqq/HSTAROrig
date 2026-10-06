"""Test with CREATE_NEW_CONSOLE flag"""
import subprocess, os, time

env = os.environ.copy()
intel_paths = [
    r'C:\Program Files (x86)\Intel\oneAPI\compiler\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\bin',
    r'C:\Program Files (x86)\Intel\oneAPI\mkl\2025.3\lib',
]
env['PATH'] = ';'.join(intel_paths) + ';' + env.get('PATH', '')

exe = r'G:\BB\hstarYLOrig\next\build_ref\hstar_ref.exe'
cwd = r'G:\BB\cube_contact\340(dizhen)'

# Test 1: CREATE_NEW_CONSOLE
print('=== Test: CREATE_NEW_CONSOLE ===')
start = time.time()
CREATE_NEW_CONSOLE = 0x00000010
result = subprocess.run(
    [exe], cwd=cwd,
    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    text=True, timeout=300, env=env,
    creationflags=CREATE_NEW_CONSOLE,
)
elapsed = time.time() - start
print(f'rc={result.returncode}, time={elapsed:.1f}s')
if result.returncode != 0:
    print(f'stderr: {result.stderr.strip()[-200:] if result.stderr else "(empty)"}')
else:
    print('SUCCESS!')

# Test 2: Use start /wait via cmd
print('\n=== Test: cmd /c start /wait ===')
start = time.time()
result = subprocess.run(
    f'cmd /c start /wait "" "{exe}"',
    cwd=cwd, shell=True, capture_output=True, text=True,
    timeout=300, env=env,
)
elapsed = time.time() - start
print(f'rc={result.returncode}, time={elapsed:.1f}s')

# Check chk file
chk = os.path.join(cwd, '1.chk')
if os.path.exists(chk):
    age = time.time() - os.path.getmtime(chk)
    with open(chk, 'r') as f:
        lines = f.readlines()
    print(f'1.chk: {len(lines)} lines, age={age:.0f}s')
