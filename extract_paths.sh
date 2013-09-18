#!/bin/bash

OUTPUT_FILE="freepath_starkexp.txt"

TMP_FILE=`mktemp`

for f in efreepath*starkexp*.txt
do
    echo -n `echo $f | sed -rn 's|efreepath_([0-9.]+)T_starkexp.txt|\1|p'` >> "$TMP_FILE"
    echo -ne "\t" >> "$TMP_FILE"
    head -n 2 "$f" | tail -n 1 | awk '{print $2}' >> "$TMP_FILE"
done

sort -g "$TMP_FILE" > "$OUTPUT_FILE"
rm "$TMP_FILE"
