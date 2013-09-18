#!/bin/bash

basedir=$1
outfile=$2
script="/home/dima/workspace/diffmc/lr_dc_window.R"

tmpfile=`mktemp`

#width="500"
#start_point="250"

width="500"
start_point="400"

for dir in "$basedir"/* ; do
    if test -d "$dir"; then

        H=`grep H "$dir"/diffmc.log | awk '{ print $4 }'`
        echo -n "$H " >> "$tmpfile"

        "$script" "$dir"/out.txt "$width" "$start_point" >> "$tmpfile"

        echo "" >> "$tmpfile"
    fi
done

#sort -g -o "$outfile" "$tmpfile"
perl -e 'print sort { $a<=>$b } <>' < "$tmpfile" > "$outfile"
rm "$tmpfile"
