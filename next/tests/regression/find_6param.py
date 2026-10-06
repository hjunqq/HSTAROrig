"""Find 6-parameter test cases and their sizes"""
import os

base_dirs = [r'G:\BB\lyr', r'G:\BB\lyc', r'G:\BB\cube_contact']
results = []

for base in base_dirs:
    for root, dirs, files in os.walk(base):
        if 'inp' in files:
            inp = os.path.join(root, 'inp')
            try:
                with open(inp, 'r') as f:
                    lines = f.readlines()
                if len(lines) >= 2:
                    params = lines[1].split('!')[0].strip().split()
                    if len(params) >= 6:
                        msh = os.path.join(root, '1.msh')
                        msh_size = os.path.getsize(msh) if os.path.exists(msh) else 0
                        chk = os.path.join(root, '1.chk')
                        chk_size = os.path.getsize(chk) if os.path.exists(chk) else 0
                        is_dizhen = 'dizhen' in root.lower()
                        results.append((msh_size, root, is_dizhen, chk_size))
            except:
                pass

results.sort()
print(f'Found {len(results)} test cases with 6 parameters:\n')
for msh_size, path, is_dizhen, chk_size in results:
    tag = ' [EARTHQUAKE]' if is_dizhen else ' [STATIC]'
    print(f'  {msh_size:>10,} bytes  {path}{tag}  chk={chk_size:>10,}')
