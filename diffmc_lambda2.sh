#!/bin/bash

H_ARRAY=("0.2e+4" "1.0e+4" "0.3e+5")
LAMBDA_ARRAY=("400.0"   "450.0"   "500.0"   "550.0"   "600.0"   "650.0"   "700.0")
TIME_ARRAY=(  "3.0e-10" "3.0e-10" "3.0e-10" "3.0e-10" "3.0e-10" "3.0e-10" "3.0e-10")


PHOTONS=100000
SCATTERINGS=10000
SEED=4000

POINTS=200
LOG_FILE="diffmc.log"

for H in "${H_ARRAY[@]}"
   do

     i=0

     for lambda in "${LAMBDA_ARRAY[@]}"
   	do

     	make clean

     	OLD_LAMBDA="const Float lambda = ;"
     	NEW_LAMBDA="const Float lambda = ${lambda}e-7;"

        OLD_H="const Float H = ;"
        NEW_H="const Float H = ${H};"

     	sed -e 's/'"$OLD_LAMBDA"'/'"$NEW_LAMBDA"'/g;s/'"$OLD_H"'/'"$NEW_H"'/g' optics.cpp.template.lambda_H > optics.cpp

     	make

     	WORK_DIR=output-$(hostname)-starkexp-${H}T-${lambda}nm_512x32768-01-e-pi_2-$(date +%Y%m%d)

     	OEPARTITION_FILE="oepartition_512x32768_${H}T_${lambda}nm_starkexp.txt"
     	EOPARTITION_FILE="eopartition_512x32768_${H}T_${lambda}nm_starkexp.txt"
     	EEPARTITION_FILE="eepartition_512x32768_${H}T_${lambda}nm_starkexp.txt"
     	OFREEPATH_FILE="ofreepath_${H}T_${lambda}nm_starkexp.txt"
     	EFREEPATH_FILE="efreepath_${H}T_${lambda}nm_starkexp.txt"
     	ECHANNELPROB_FILE="echannel_${H}T_${lambda}nm_starkexp.txt"

     	MAXTIME=${TIME_ARRAY[$i]}

     	mkdir -p $WORK_DIR

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
          2> $WORK_DIR/$LOG_FILE

   	 echo "point H=${H}, lambda=${lambda}nm done"
    	((i++))
     done
  done
