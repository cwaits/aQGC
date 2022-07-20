#!/bin/bash -x

dir=$1
touch "$dir"/limitsTable
for cut_combo in "$dir"*; do
    if [[ "$cut_combo" != "$dir""limitsTable" ]];
    then
        echo "$cut_combo"
        #cat "$cut_combo"/log_limits
        limits="$(grep  95% $cut_combo/log_limits | tail -1)"
        cut_dict="$(grep cut_dict $cut_combo/cut_combo.py)"
        cut="$(basename -- $cut_combo)"
        table="$dir""limitsTable"
        entry="$cut_dict: $limits"
        echo "$cut:" >> "$table"
        echo "$entry" >> "$table"
        #echo " " >> "$table"
    fi
done
python2.7 sortLimits.py "$dir""limitsTable"
