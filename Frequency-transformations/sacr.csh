#!/bin/csh
ls *.sac | wc -l > log
ls *.sac >> log
./sacr.e
rm -rf log
