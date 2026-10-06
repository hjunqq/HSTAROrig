import sys
with open(r'G:\BB\orig\7static\S1\inp', 'rb') as f:
    data = f.read()
print(f'Size: {len(data)} bytes')
print(f'Content: {repr(data[:200])}')
