#!/bin/bash

basedir=$1
outfile=$2

tmpfile=`mktemp`
data=

for file in "$basedir"/* ; do

    if ["$data"] then
        data="$tmpfile"
    else
        data=file
    fi

    if test -f "$dir"; then
        join file data > "$tmpfile"
    fi
done

mv "$tmpfile" "$outfile"
