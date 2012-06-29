#!/bin/bash

# plotR4RNA.sh

WIDTH=600
HEIGHT=300

for FOLD in *.fold
do
    echo $FOLD
    PNG=${FOLD%fold}png
    echo "library(R4RNA)" > R.scr
    echo "rna <- readVienna(\"$FOLD\")" >> R.scr
    echo "png(\"$PNG\", width=$WIDTH, height=$HEIGHT)" >> R.scr
    echo "plotHelix(rna)" >> R.scr
    R --vanilla --quiet < R.scr
done

rm R.scr
