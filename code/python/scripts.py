#!/usr/bin/env python 

import argparse
import os

parser = argparse.ArgumentParser(description='Process inversions, repeats and recombination maps and divides genome into windows.')
parser.add_argument("--inversions", type=str, metavar = "FILE", help="inversions csv")
parser.add_argument("--repeats", type=str,metavar="FILE", help="file with repeats")

args = (l,b)
(l,b) = [ value for _, value in parser.parse_args()._get_kwargs() ]
print(args)