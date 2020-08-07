#!/usr/bin/env python
import ph4
import sys,os

usage = """
usage: make_one_ph4.py <table.txt> <min elements to match> <title>

e.g.

make_one_ph4.py hivp_nmr.txt 1 "This is my table" > hivp_nmr.ph4

(don't forget the quotes around the title)
"""
if len(sys.argv) != 4:
    sys.exit(usage)
print ph4.get_ph4_from_table(file(sys.argv[1]).read(),sys.argv[2],sys.argv[3])
