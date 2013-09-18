#!/bin/bash

WORK_DIR=output-$(hostname)-starkexp-1.0T-4096x4096-01-e-pi_2-$(date +%Y%m%d)
#OEPARTITION_FILE="oepartition01.txt"
#EOPARTITION_FILE="eopartition01.txt"
#EEPARTITION_FILE="eepartition01.txt"
#OFREEPATH_FILE="ofreepath04.txt"
#EFREEPATH_FILE="efreepath04.txt"
#ECHANNELPROB_FILE="echannel04.txt"

OEPARTITION_FILE="oepartition_1.0T_starkexp.txt"
EOPARTITION_FILE="eopartition_1.0T_starkexp.txt"
EEPARTITION_FILE="eepartition_1.0T_starkexp.txt"
OFREEPATH_FILE="ofreepath_1.0T_4096x4096_starkexp.txt"
EFREEPATH_FILE="efreepath_1.0T_4096x4096_starkexp.txt"
ECHANNELPROB_FILE="echannel_1.0T_4096x4096_starkexp.txt"


#PHOTONS=5000
PHOTONS=10000
SCATTERINGS=10000
#SCATTERINGS=1000
SEED=2000
MAXTIME='1e-10'
#MAXTIME='5e-9'
#MAXTIME='1'
#POINTS=500
H='0.2'
POINTS=1000
LOG_FILE="diffmc.log"

mkdir -p $WORK_DIR

#gdb -ex run --args  \
./diffmc \
        --workdir $WORK_DIR \
        --saveoepartition $OEPARTITION_FILE \
        --saveeopartition $EOPARTITION_FILE \
        --saveeepartition $EEPARTITION_FILE \
        --saveofreepath $OFREEPATH_FILE \
        --saveefreepath $EFREEPATH_FILE \
        --saveechannelprob $ECHANNELPROB_FILE \
        --photons $PHOTONS \
        --scatterings $SCATTERINGS \
        --maxtime $MAXTIME \
        --points $POINTS\
        --seed $SEED \
        2> $WORK_DIR/$LOG_FILE &

tail -f $WORK_DIR/$LOG_FILE
