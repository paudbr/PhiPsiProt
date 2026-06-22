#!/usr/bin/env python3
import sys
seq = ''
with open(sys.argv[1]) as f:
    for line in f:
        if not line.startswith('>'):
            seq += line.strip()
print(seq, end='')