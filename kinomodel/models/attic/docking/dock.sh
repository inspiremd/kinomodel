#!/bin/bash
#BSUB -J shapefit
#BSUB -n 1
#BSUB -R span[ptile=1]
#BSUB -R rusage[mem=20]
#BSUB -W 20:00
#BSUB -o %J.out
#BSUB -eo %J.err
 
cd $LS_SUBCWD
python shapefit.py
