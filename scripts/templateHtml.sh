#!/bin/sh

#define parameters which are passed in.
RUNNUMBER=$1

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
    <h2>Last run number is $RUNNUMBER</h2>
    <img src="plots/Run${RUNNUMBER}_hitMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
							margin-bottom:
							0.5em;">

    <img src="plots/Run${RUNNUMBER}_coincidenceMap.png" style="float: left;
							width: 47%;
							margin-right:
							1%;
							margin-bottom:
							0.5em;">


  </div>


</body>
</html>
EOF