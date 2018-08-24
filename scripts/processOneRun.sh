#!/bin/bash

source /home/lindac/DUNE/TOF/DsTOFanalysis/env.sh

cd /home/lindac/DUNE/TOF/DsTOFanalysis/build/

export DATADIR=/unix/dune/hptpctof/
export htmlPageDir=/home/lindac/public_html/dstof/

runNumber=$1;

COUNTER=0

chmod 777 $DATADIR/run$runNumber
chmod 777 $DATADIR/run$runNumber/*
hitmap=$DATADIR/run$runNumber/Run${runNumber}_histos.root
niceTree=$DATADIR/run${runNumber}/DsTOFtreeRun${runNumber}_tdc2.root


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
