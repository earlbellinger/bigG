#!/bin/bash
#### Script for dispatching an ASTEC simulation 
#### Author: Earl Bellinger ( bellinger@phys.au.dk ) 
#### Stellar Astrophysics Centre Aarhus 

# run evolution code 
source astec.sh 

MDL_NUMS=$(tail -n +5 LOGS/history.data | awk '{ print $1 }')

for MDL_NUM in $MDL_NUMS; do
    # obtain the structure 
    csh run-evol.mass.stg d.$MT $ZT $X $ALPHAT $OV $CASE $MDL_NUM 
    
    # redistribute the mesh 
    csh run.redistrb prxq4 d.$MT.Z$ZT.X$X.a$ALPHAT.o$OV.$CASE $MDL_NUM 
    
    # calculate adiabatic pulsation frequencies 
    csh run.adipls d.$MT.Z$ZT.X$X.a$ALPHAT.o$OV.$CASE $MDL_NUM prxq4 
    
    # process frequencies 
    FREQ_FILE="osc/fobs.d.$MT.Z$ZT.X$X.a$ALPHAT.o$OV.$CASE.$MDL_NUM.prxq4"
    #cat $FREQ_FILE 
    cp $FREQ_FILE LOGS/profile$MDL_NUM-freqs.dat
    
    FGONG_FILE="gong/fgong.d.$MT.Z$ZT.X$X.a$ALPHAT.o$OV.$CASE.$MDL_NUM"
    FGONG=LOGS/profile$MDL_NUM.FGONG
    cp $FGONG_FILE $FGONG
    python3 /home/earl/asteroseismology/scripts/fgong2ascii.py -i $FGONG
done

cd - 
Rscript summarize-track.R "$DIRNAME"
if [ $REMOVE -eq 1 ]; then
    rm -rf "$DIRNAME"
fi

# fin. 
