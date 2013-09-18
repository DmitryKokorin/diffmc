#!/bin/bash

basedir=$1
outfile=$2
script="/home/dima/workspace/diffmc/ml.R"

tmpfile=`mktemp`

for dir in "$basedir"/* ; do
    if test -d "$dir"; then
        H=`grep H "$dir"/diffmc.log | awk '{ print $4 }'`
        maxtime=`grep maxtime "$dir"/diffmc.log | awk '{ print $4 }'`
        "$script" "$dir"/pos.txt "$H" "$maxtime" >> "$tmpfile"
    fi
done

sort -g -o "$outfile" "$tmpfile"
rm tmpfile
