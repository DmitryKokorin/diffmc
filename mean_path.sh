#!/bin/bash

for f in efreepath_*T_starkexp_3_ee_only_equal_Kll.txt ; do
    ./mean_path.py $f "mean_$f"
    #echo $f
done
