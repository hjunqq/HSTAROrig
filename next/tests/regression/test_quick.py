"""Quick test: does hstar_ref.exe with /Qinit:zero work?"""
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

print(f'Running {os.path.basename(exe)}...')
start = time.time()
result = subprocess.run(
    [exe], cwd=cwd, capture_output=True, text=True, timeout=300, env=env,
)
elapsed = time.time() - start
print(f'rc={result.returncode}, time={elapsed:.1f}s')
if result.returncode != 0:
    for line in result.stderr.strip().splitlines()[-5:]:
        print(f'  {line}')
else:
    print('SUCCESS!')
    stdout_lines = result.stdout.strip().splitlines()
    print(f'stdout: {len(stdout_lines)} lines')
    for line in stdout_lines[-5:]:
        print(f'  {line}')
