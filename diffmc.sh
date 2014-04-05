#!/bin/bash

#DATA_NAME=stark_lubensky_pre_1997
#DATA_NAME=van_tiggelen_stark_rmp_2000_1
#DATA_NAME=van_tiggelen_stark_rmp_2000_2
#DATA_NAME=van_tiggelen_mclc_1997_paa
#DATA_NAME=van_tiggelen_mclc_1997_mbba_prec0005_minwidth_1e-6
#DATA_NAME=heiderich_1997_1
#DATA_NAME=heiderich_1997_2
#DATA_NAME=kao_1997_v2
#DATA_NAME=wiersma_1999_v2
#DATA_NAME="figures_start_h=9.0T"
DATA_NAME="figures_v3_h=0.2T_angle"

WORK_DIR=output-$(hostname)-${DATA_NAME}-e-pi_3-$(date +%Y%m%d)

OEPARTITION_FILE="oepartition_${DATA_NAME}.txt"
EOPARTITION_FILE="eopartition_${DATA_NAME}.txt"
EEPARTITION_FILE="eepartition_${DATA_NAME}.txt"
OFREEPATH_FILE="ofreepath_${DATA_NAME}.txt"
EFREEPATH_FILE="efreepath_${DATA_NAME}.txt"
ECHANNELPROB_FILE="echannel_${DATA_NAME}.txt"


PHOTONS=5000000
#PHOTONS=1000000
#PHOTONS=100000
SCATTERINGS=1000000
SEED=43
#MAXTIME='2e-9'
MAXTIME='2e-10'
#MAXTIME='5e-9'
#MAXTIME='1'
#POINTS=500
POINTS=1000
LOG_FILE="diffmc.log"

mkdir -p $WORK_DIR

#gdb -ex run --args  \
./diffmc \
        --workdir $WORK_DIR \
        --loadoeprofile $OEPARTITION_FILE \
        --loadeoprofile $EOPARTITION_FILE \
        --loadeeprofile $EEPARTITION_FILE \
        --loadofreepath $OFREEPATH_FILE \
        --loadefreepath $EFREEPATH_FILE \
        --loadeeprobability $ECHANNELPROB_FILE \
        --photons $PHOTONS \
        --scatterings $SCATTERINGS \
        --maxtime $MAXTIME \
        --points $POINTS\
        --seed $SEED \
        2> $WORK_DIR/$LOG_FILE &

tail -f $WORK_DIR/$LOG_FILE
