#!/bin/sh
# properties = {properties}
hostname=`hostname`
echo "Running on $hostname"
startTime=`date`
echo "Start time: $startTime"
{exec_job}
endTime=`date`
echo "End time: $endTime"
