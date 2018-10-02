#/bin/bash

find . -name "*.smi" | while read f
do
  awk '{print $2}' $f > ${f%sim}.sim1
done
