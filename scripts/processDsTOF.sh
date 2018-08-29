#!/bin/bash

source /home/lindac/DUNE/TOF/DsTOFanalysis/env.sh

cd /home/lindac/DUNE/TOF/DsTOFanalysis/build/

export DATADIR=/unix/dune/hptpctof/
export htmlPageDir=/home/lindac/public_html/dstof/

LASTRUN=`ls -t /unix/dune/hptpctof/run*/raw* | head -n 1`
LASTRUN=${LASTRUN#*run}
LASTRUN=${LASTRUN%/raw*}
echo $LASTRUN

COUNTER=0
for dir in `ls $DATADIR`; 
do
    echo $dir

    if [[ $dir =~ [0-9] ]]; then
	runNumber=$(echo $dir | sed 's/[^0-9]*//g')
        echo 'Run number is' $runNumber
    else
	continue
    fi

    if [ "$runNumber" -lt 20 ]; then
	continue
    fi

    if [ "$runNumber" == 128 ]; then
	continue
    fi

    if [ "$runNumber" == 129 ]; then
	continue
    fi

    if [ "$runNumber" == 130 ]; then
	continue
    fi

    if [ "$runNumber" == 145 ]; then
	continue
    fi

    if [ "$runNumber" == 724 ]; then
	continue
    fi

    if [ "$runNumber" == 758 ]; then
	continue
    fi

    chmod 777 $DATADIR/run$runNumber
    chmod 777 $DATADIR/run$runNumber/*

    if [ "$runNumber" == $LASTRUN ]; then
	continue
    fi
    
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

    cp $DATADIR/run$runNumber/Run${runNumber}_hitMap.png                 $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_hitTimeMap.png             $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_hitBeamSpillMap.png        $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceMap.png         $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_barEfficiency.png          $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceInSpill.png     $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceInSpillBar.png  $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_hitsInSpillBar.png         $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_Efficiency.png             $htmlPageDir/plots/
    
    
    
done

./makeHVplots $htmlPageDir/

echo "Last rootified run number is " $[$LASTRUN -1 ]

/home/lindac/DUNE/TOF/DsTOFanalysis/scripts/templateHtml.sh $[$LASTRUN -1 ] > $htmlPageDir/dstofSummary.html

