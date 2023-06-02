#!/bin/sh
# Give the job a name
$ -N JOB_DOL
# set the shell
$ -S /bin/sh
# set working directory on all host to directory where the job was started
$ -cwd
# send all ERROR messages to this file
$ -e errors.txt
# Change the email address to YOUR email, and you will be emailed when the job has finished.
$ -m e
$ -M bankheak@oregonstate.edu
# Ask for 1 core, as R can only use 1 core for processing
$ -pe orte 1
# Load the R Module
module load R
# Commands to run job
R --slave < breakdown_test.r > output.${SGE_TASK_ID} --args $mystart $myend
