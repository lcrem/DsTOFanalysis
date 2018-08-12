#!/bin/bash

source /home/lindac/DUNE/TOF/DsTOFanalysis/env.sh

cd /home/lindac/DUNE/TOF/DsTOFanalysis/build/

export DATADIR=/unix/dune/hptpctof/
export htmlPageDir=/home/lindac/public_html/dstof/

for dir in `ls $DATADIR`; 
do
    echo $dir
    if [[ $dir =~ [0-9] ]]; then
	runNumber=$(echo $dir | sed 's/[^0-9]*//g')
        echo 'Run number is' $runNumber
    else
	continue
    fi

    if [ "$runNumber" -lt 2 ]; then
	continue
    fi
    

    hitmap=$DATADIR/run$runNumber/Run${runNumber}_histos.root
    niceTree=$DATADIR/run${runNumber}/DsTOFtreeRun${runNumber}_tdc2.root

    if [ -f $hitmap ]; then
	echo "Found Everything!"
	echo "Cleaning up parsed files"
	if  ls $DATADIR/run${runNumber}/parsed* 1> /dev/null 2>&1; then
	    echo "Parsed files found!"
	    for parsed in `$DATADIR/run${runNumber}/parsed*`; do
		echo "Removing parsed file " $parsed
		rm $parsed
	    do
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
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceMap.png  $htmlPageDir/plots/
    
    
    
done

echo "Last rootified run number is " $runNumber

/home/lindac/DUNE/TOF/DsTOFanalysis/scripts/templateHtml.sh $runNumber > $htmlPageDir/dstofSummary.html

