#!/bin/bash

source /home/lindac/DUNE/TOF/DsTOFanalysis/env.sh

cd /home/lindac/DUNE/TOF/DsTOFanalysis/build/

export DATADIR=/unix/dune/hptpctof/
export htmlPageDir=/home/lindac/public_html/dstof/

TOTAL=0
for dir in `ls $DATADIR`; 
do
 TOTAL=$[$TOTAL +1]
done

echo "Total number of dirs is " $TOTAL


COUNTER=0
for dir in `ls -rt $DATADIR`; 
do
    echo $dir
    COUNTER=$[$COUNTER +1]
    if [ $COUNTER -ge $TOTAL ]; then
	echo "Skipping the last dir"
	continue
    fi


    if [[ $dir =~ [0-9] ]]; then
	runNumber=$(echo $dir | sed 's/[^0-9]*//g')
        echo 'Run number is' $runNumber
    else
	continue
    fi

    if [ "$runNumber" -lt 20 ]; then
	continue
    fi
    

    chmod 777 $DATADIR/run$runNumber
    chmod 777 $DATADIR/run$runNumber/*
    hitmap=$DATADIR/run$runNumber/Run${runNumber}_histos.root
    niceTree=$DATADIR/run${runNumber}/DsTOFtreeRun${runNumber}_tdc2.root

    if [ -f $hitmap ]; then
	echo "Found Everything!"
	echo "Cleaning up parsed files"
	if  ls $DATADIR/run${runNumber}/parsed* 1> /dev/null 2>&1; then
	    echo "Parsed files found!"
	    for parsed in `ls $DATADIR/run${runNumber}/parsed*`; do
		ls $parsed
		echo "Removing parsed file " $parsed
		rm $parsed
	    done
	    
	fi
	continue
    fi
    

    if [ -f $niceTree ]; then
	echo "Nice trees exist!"    
    else 
	echo "Producing nice tree .."

	echo "Check if parsed files exist:"
	if  ls $DATADIR/run${runNumber}/parsed* 1> /dev/null 2>&1; then
	    echo "Parsed files found!"
	else
	    echo "Parsed files not found"
	    cd /home/lindac/DUNE/TOF/hptpcTofVmeReadout/
	    ./decodeRun.sh $runNumber
	    cd -
	fi

	./makeNiceTree $runNumber
	
    fi
    
    echo "Producing map tree!"
    ./makeHitMap $runNumber
    cp $DATADIR/run$runNumber/Run${runNumber}_hitMap.png          $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_hitTimeMap.png      $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_hitBeamSpillMap.png $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceMap.png  $htmlPageDir/plots/
    
    
    
done

echo "Last rootified run number is " $runNumber

/home/lindac/DUNE/TOF/DsTOFanalysis/scripts/templateHtml.sh $runNumber > $htmlPageDir/dstofSummary.html

