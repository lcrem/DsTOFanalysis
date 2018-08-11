#!/bin/bash

source /home/lindac/DUNE/TOF/DsTOFanalysis/env.sh

cd /home/lindac/DUNE/TOF/DsTOFanalysis/build/

export DATADIR=/unix/dune/hptpctof/
htmlPageDir=/home/lindac/public_html/dstof/

for dir in `ls $DATADIR`; 
do
    echo $dir
    if [[ $dir =~ [0-9] ]]; then
	runNumber=$(echo $dir | sed 's/[^0-9]*//g')
        echo 'Run number is' $runNumber
    else
	continue
    fi

    if [ "$runNumber" -lt 25 ]; then
	continue
    fi
    
    
    echo "Check if parsed files exist:"

    if ls $dir/parsed*txt 1> /dev/null 2>&1; then
	echo "Parsed files not found"
    else
	echo "Parsed files found"
	niceTree=$dir/DsTOFtree_tdc1.root
	if [ ! -f $niceTree ]; then
	    echo "Producing nice tree"
	    ./makeNiceTree $runNumber
	    else
	    echo "Nice trees exist!"
	fi
	./makeHitMap $runNumber
	
	cp $DATADIR/run$runNumber/Run${runNumber}_hitMap.png          $htmlPageDir/plots/
	cp $DATADIR/run$runNumber/Run${runNumber}_coincidenceMap.png  $htmlPageDir/plots/
	
    fi
    
done

echo "Last rootified run number is " $runNumber

/home/lindac/DUNE/TOF/DsTOFanalysis/scripts/templateHtml.sh $runNumber > $htmlPageDir/dstofSummary.html
