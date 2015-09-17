#!/bin/bash

DATE=`date +%Y%m%d`
if [[ `uname` == 'Linux' ]]; then
    CSVDIR="/home/$USER/testSuiteTiming"
else
    CSVDIR="/Users/$USER/testSuiteTiming"
fi

if [ ! -d $CSVDIR ]; then
    mkdir -p $CSVDIR
fi

for file in $1 $2; do
    time=`grep cputime $file | cut -d: -f3 | sed -e 's/^[[:space:]]*//' | sed -e 's/[[:space:]]*$//'`
    csvfile=`echo $file | sed 's/\//-/g' | sed 's/\.trs$/\.csv/g'`
    if [ ! -f $CSVDIR/$csvfile ]; then
        name=`echo $file | sed 's/\//-/g' | sed 's/\.m.trs$//g' | sed 's/\.o.trs$//g'`
        echo "DATE,$name" > $CSVDIR/$csvfile
    fi
    echo $DATE,$time >> $CSVDIR/$csvfile
done
