import sys
with sys.stdin as inf:
    for line in inf:
        if not line.startswith('#'):
            if not any(x in line.split()[0] for x in ['X','Y','M']): sys.stdout.write(line)
        else: sys.stdout.write(line)
