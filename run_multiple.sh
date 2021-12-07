#!/bin/sh

for fSize in 1 2 4 8 16 32
do
  python3 Parallel.py -p=8 -l=20 -g=3 -d='multiple/synthesis_'$fSize'_6219' -e=0 -f=1
done