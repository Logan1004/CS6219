#!/bin/sh

#This program returns the
#contents of my Home folder

for i in 1 2 4 8 16 32
do
  echo 'synthesis_'$i'_6219'
  python3 Utils.py  'multiple/synthesis_'$i'_6219/Primer/primer.txt' 'multiple/synthesis_'$i'_6219/Primer/primer_gram'
done