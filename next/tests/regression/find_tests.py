"""Find all inp files with content and check analysis types"""
import os
from pathlib import Path

base_dirs = [r'G:\BB\orig', r'G:\BB\lyr', r'G:\BB\static', r'G:\BB\lyc']

for base in base_dirs:
    for root, dirs, files in os.walk(base):
        if 'inp' in files:
            inp_path = os.path.join(root, 'inp')
            size = os.path.getsize(inp_path)
            if size == 0:
                continue
            try:
                with open(inp_path, 'r', errors='replace') as f:
                    lines = f.readlines()
                print(f'\n{root}')
                print(f'  Size: {size} bytes, {len(lines)} lines')
                for i, line in enumerate(lines[:6]):
                    print(f'  [{i+1}] {line.rstrip()}')
                # Check for .msh file
                msh_files = [f for f in files if f.endswith('.msh') and not f.endswith('.flavia.msh')]
                print(f'  .msh files: {msh_files[:3]}')
            except Exception as e:
                print(f'  ERROR: {e}')
