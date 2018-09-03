#!/bin/sh

#define parameters which are passed in.
RUNNUMBER=$1
LASTRUN=$[$RUNNUMBER -5]

#define the template.
cat  <<EOF
<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">

  <title>High Pressure TPC Downstream TOF</title>
  <meta http-equiv="refresh" content="300" >
  <meta name="description" content="High Pressure TPC Downstream TOF">
  <meta name="author" content="LindaC">

</head>

<body>

  <div>
    <h1>High Pressure TPC Downstream TOF</h1>
    <h2>Last rootified run is $RUNNUMBER. Data shouldn't be older than 1.5 hours!</h2>
    <a href="hvstatus.html">Click here to check HV status</a>
 and <a href="hvplots.html"> here to check HV plots</a> (updates every 15 mins)
</div>
EOF

for ((run=$RUNNUMBER;run>LASTRUN;run-=1)) ; do 
    cat  <<EOF
<div>
    <h2>Plots for run $run</h2>
    <img src="plots/Run${run}_hitMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
							margin-bottom:
							0.5em;">
    <img src="plots/Run${run}_hitTimeMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
							margin-bottom:
							0.5em;">

    <img src="plots/Run${run}_coincidenceMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
							margin-bottom:
							0.5em;">

    <img src="plots/Run${run}_hitBeamSpillMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_barEfficiency.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_coincidenceInSpill.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_coincidenceInSpillBar.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_hitsInSpillBar.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_Efficiency.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">
    <img src="plots/Run${run}_ToF.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
  							margin-bottom:
						0.5em;">

  </div>
EOF
done

cat  <<EOF
  </div>
</body>
</html>
EOF