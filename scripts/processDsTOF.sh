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

    if [ "$runNumber" -lt 26 ]; then
	continue
    fi
    
    
    echo "Check if parsed files exist:"
    if  ls $DATADIR/run${runNumber}/parsed* 1> /dev/null 2>&1; then
	echo "Parsed files found!"
    else
	echo "Parsed files not found"
	cd /home/lindac/DUNE/TOF/hptpcTofVmeReadout/
	./decodeRun.sh $runNumber
	cd -
    fi

    niceTree=$dir/DsTOFtree_tdc1.root
    if [ -f $niceTree ]; then
	echo "Nice trees exist!"    
    else 
	echo "Producing nice tree"
	./makeNiceTree $runNumber
    fi

    ./makeHitMap $runNumber
    cp $DATADIR/run$runNumber/Run${runNumber}_hitMap.png          $htmlPageDir/plots/
    cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceMap.png  $htmlPageDir/plots/
    
done

echo "Last rootified run number is " $runNumber

/home/lindac/DUNE/TOF/DsTOFanalysis/scripts/templateHtml.sh $runNumber > $htmlPageDir/dstofSummary.html
