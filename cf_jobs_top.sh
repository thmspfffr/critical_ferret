3#!/bin/sh

# Submits n jobs to the torque queing system

for i in {1..20}
do
  echo 'Start Job' $i
  qsub cf_jobs.sh
  sleep 10
done
