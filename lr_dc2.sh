#!/bin/bash

basedir=$1
outfile=$2
script="/home/dima/workspace/diffmc/lr_dc2.R"

tmpfile=`mktemp`

width="1000"
start_array=["1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000"]

for dir in "$basedir"/* ; do
    if test -d "$dir"; then

        H=`grep H "$dir"/diffmc.log | awk '{ print $4 }'`
        echo "$H" >> "$tmpfile"

        for start_value in "$start_array[@]"
        do
            echo " " >> "$tmpfile"
            "$script" "$dir"/out.txt "$width" "$start_value" >> "$tmpfile"
        done

        echo "\n" >> "$tmpfile
        "
    fi
done

sort -g -o "$outfile" "$tmpfile"
rm "$tmpfile"
