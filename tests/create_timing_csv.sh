#!/bin/bash

DATE=`date +%Y%m%d`
CSVDIR="/home/dynbot/testSuiteTiming"
if [ ! -d $CSVDIR ]; then
    mkdir -p $CSVDIR
fi

for file in $1 $2; do
    time=`grep cputime $file | cut -d: -f3 | sed -e 's/^[[:space:]]*//' | sed -e 's/[[:space:]]*$//'`
    csvfile=`echo $file | sed 's/\//-/g' | sed 's/\.trs$/\.csv/g'`
    if [ ! -f $CSVDIR/$csvfile ]; then
        `touch $CSVDIR/$csvfile`
    fi
    echo $DATE,$time >> $CSVDIR/$csvfile
done
