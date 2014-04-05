#!/bin/bash

#H_ARRAY=("0.2"        "0.3"      "0.4"      "0.5"      "0.6"      "0.7"      "0.8"      "0.9"      "1.0"      "2.0"      "3.0"      "6.0"     "9.0"      "12.0"    "15.0"    "18.0"     "21.0"    "24.0"    "27.0"     "30.0"    "33.0"    "36.0")
#TIME_ARRAY=("2.44e-8" "2.44e-8"  "2.54e-8"  "2.66e-8"  "2.66e-8"  "2.70e-8"  "2.80e-8"  "3.0e-8"   "3.0e-8"   "3.18e-8"  "3.38e-8"  "4.0e-8"  "4.76e-8"  "5.6e-8"  "6.0e-8"  "6.66e-8"  "7.6e-8"  "8.0e-8"  "8.68e-8"  "9.6e-8"  "1.04e-7" "1.102e-7")
#TIME_ARRAY=("2.44e-10" "2.44e-10" "2.54e-10" "2.66e-10" "2.66e-10" "2.70e-10" "2.80e-10" "3.0e-10" "3.0e-10" "3.18e-10" "3.38e-10" "4.0e-10" "4.76e-10" "5.6e-10" "6.0e-10" "6.66e-10" "7.6e-10" "8.0e-10" "8.68e-10" "9.6e-10" "1.04e-9" "1.102e-9")

#H_ARRAY=("4.0"       "5.0"    "7.0"    "8.0"    "10.0"   "11.0"   "13.0"    "14.0"    "16.0"   "17.0"   "19.0"   "20.0"   "22.0"   "23.0"   "25.0"   "26.0"   "28.0"   "29.0"   "31.0"   "32.0"   "34.0"    "35.0")
#TIME_ARRAY=("3.5e-9" "3.7e-9" "4.3e-9" "4.6e-9" "5.0e-9" "5.3e-9" "5.75e-9" "5.85e-9" "6.2e-9" "6.4e-9" "7.0e-9" "7.3e-9" "7.8e-9" "7.9e-9" "8.2e-9" "8.4e-9" "9.0e-9" "9.3e-9" "9.8e-9" "1.0e-8" "1.06e-8" "1.08e-8")

#H_ARRAY=("0.2"        "0.3"      "0.4"      "0.5"      "0.6"      "0.7"      "0.8"      "0.9"      "1.0"      "2.0"      "3.0"      "6.0"     "9.0"      "12.0"    "15.0"    "18.0"     "21.0"    "24.0"    "27.0"     "30.0"    "33.0"    "36.0"     "4.0"       "5.0"    "7.0"    "8.0"    "10.0"   "11.0"   "13.0"    "14.0"    "16.0"   "17.0"   "19.0"   "20.0"   "22.0"   "23.0"   "25.0"   "26.0"   "28.0"   "29.0"   "31.0"   "32.0"   "34.0"    "35.0")
#TIME_ARRAY=("2.44e-9" "2.44e-9"  "2.54e-9"  "2.66e-9"  "2.66e-9"  "2.70e-9"  "2.80e-9"  "3.0e-9"   "3.0e-9"   "3.18e-9"  "3.38e-9"  "4.0e-9"  "4.76e-9"  "5.6e-9"  "6.0e-9"  "6.66e-9"  "7.6e-9"  "8.0e-9"  "8.68e-9"  "9.6e-9"  "1.04e-8" "1.102e-8" "3.5e-9" "3.7e-9" "4.3e-9" "4.6e-9" "5.0e-9" "5.3e-9" "5.75e-9" "5.85e-9" "6.2e-9" "6.4e-9" "7.0e-9" "7.3e-9" "7.8e-9" "7.9e-9" "8.2e-9" "8.4e-9" "9.0e-9" "9.3e-9" "9.8e-9" "1.0e-8" "1.06e-8" "1.08e-8")

#H_ARRAY=("14.0")
#TIME_ARRAY=("5.85e-9")

H_ARRAY=("0.2"         "0.5"     "1.0"      "2.0"      "3.0"      "6.0"     "9.0"      "12.0"    "15.0"    "18.0"     "21.0"    "24.0"    "27.0"     "30.0"    "33.0"    "36.0"    "4.0"       "5.0"    "7.0"    "8.0"    "10.0"   "11.0"   "13.0"    "14.0"    "16.0"   "17.0"   "19.0"   "20.0"   "22.0"   "23.0"   "25.0"   "26.0"   "28.0"   "29.0"   "31.0"   "32.0"   "34.0"    "35.0")
TIME_ARRAY=("2.44e-9"  "2.66e-9" "3.0e-9"   "3.18e-9"  "3.38e-9"  "4.0e-9"  "4.76e-9"  "5.6e-9"  "6.0e-9"  "6.66e-9"  "7.6e-9"  "8.0e-9"  "8.68e-9"  "9.6e-9"  "1.04e-8" "1.102e-8" "3.5e-9" "3.7e-9" "4.3e-9" "4.6e-9" "5.0e-9" "5.3e-9" "5.75e-9" "5.85e-9" "6.2e-9" "6.4e-9" "7.0e-9" "7.3e-9" "7.8e-9" "7.9e-9" "8.2e-9" "8.4e-9" "9.0e-9" "9.3e-9" "9.8e-9" "1.0e-8" "1.06e-8" "1.08e-8")

#H_ARRAY=("0.2")
#TIME_ARRAY=("2.44e-9")
#H_ARRAY=("15.0")
#TIME_ARRAY=("6.0e-9")


POINTS=1000
#POINTS=20000
PHOTONS=1000000
#PHOTONS=5000000
SCATTERINGS=1000000
#SEED=42
#SEED=43
SEED=44

LOG_FILE="diffmc.log"

i=0

for h in "${H_ARRAY[@]}"
   do

     make clean

     OLD="const Float H = ;"
     NEW="const Float H = ${h}e+4;"

     sed -e 's/'"$OLD"'/'"$NEW"'/g' optics.cpp.template.H.test > optics.cpp
     make

#     WORK_DIR=output-$(hostname)-starkexp-${h}T-32768x32768-01-e-pi_2-$(date +%Y%m%d)

#     OEPROFILE_FILE="oeprofile_${h}T_512x32768_starkexp.txt"
#     EOPROFILE_FILE="eoprofile_${h}T_512x32768_starkexp.txt"
#     EEPROFILE_FILE="eeprofile_${h}T_512x32768_starkexp.txt"
#     OFREEPATH_FILE="ofreepath_${h}T_starkexp.txt"
#     EFREEPATH_FILE="efreepath_${h}T_starkexp.txt"
##     ECHANNELPROB_FILE="echannel_${h}T_starkexp.txt"
#     ECHANNELPROB_FILE="echannel_ee_only.txt"

     WORK_DIR=output-$(hostname)-test-${h}T-0.01-2-e-pi_2-$(date +%Y%m%d)

     OEPROFILE_FILE="oeprofile_dummy.txt"
     EOPROFILE_FILE="eoprofile_dummy.txt"
     EEPROFILE_FILE="eeprofile_${h}T_0.01_test_2.txt"
     OFREEPATH_FILE="ofreepath_dummy.txt"
     EFREEPATH_FILE="efreepath_${h}T_test_2.txt"
     ECHANNELPROB_FILE="echannel_ee_only.txt"

     MAXTIME=${TIME_ARRAY[$i]}

     mkdir -p $WORK_DIR

     if [ -f ${OEPROFILE_FILE} ]
     then
        OEPROFILE_OPTION="--loadoeprofile ${OEPROFILE_FILE}"
     else
        OEPROFILE_OPTION="--saveoeprofile ${OEPROFILE_FILE}"
     fi

     if [ -f ${EOPROFILE_FILE} ]
     then
        EOPROFILE_OPTION="--loadeoprofile ${EOPROFILE_FILE}"
     else
        EOPROFILE_OPTION="--saveeoprofile ${EOPROFILE_FILE}"
     fi

     if [ -f ${EEPROFILE_FILE} ]
     then
        EEPROFILE_OPTION="--loadeeprofile ${EEPROFILE_FILE}"
     else
        EEPROFILE_OPTION="--saveeeprofile ${EEPROFILE_FILE}"
     fi

     if [ -f ${OFREEPATH_FILE} ]
     then
        OFREEPATH_OPTION="--loadofreepath ${OFREEPATH_FILE}"
     else
        OFREEPATH_OPTION="--saveofreepath ${OFREEPATH_FILE}"
     fi

     if [ -f ${EFREEPATH_FILE} ]
     then
        EFREEPATH_OPTION="--loadefreepath ${EFREEPATH_FILE}"
     else
        EFREEPATH_OPTION="--saveefreepath ${EFREEPATH_FILE}"
     fi

     if [ -f ${ECHANNELPROB_FILE} ]
     then
        ECHANNELPROB_OPTION="--loadeeprobability ${ECHANNELPROB_FILE}"
     else
        ECHANNELPROB_OPTION="--saveeeprobability ${ECHANNELPROB_FILE}"
     fi

#    valgrind \
#     gdb -ex run --args \
     ./diffmc \
          --workdir $WORK_DIR \
          $OEPROFILE_OPTION \
          $EOPROFILE_OPTION \
          $EEPROFILE_OPTION \
          $OFREEPATH_OPTION \
          $EFREEPATH_OPTION \
          $ECHANNELPROB_OPTION \
          --photons $PHOTONS \
          --scatterings $SCATTERINGS \
	  --maxtime $MAXTIME \
          --points $POINTS\
          --seed $SEED \
          2> $WORK_DIR/$LOG_FILE

    echo "point H=${h}T done"
    ((i++))
   done
