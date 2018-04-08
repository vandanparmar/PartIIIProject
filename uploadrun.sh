#!/bin/bash
scp *.py $1:~
scp *.json $1:~
ssh -f $1 "python3 batch.py>>log.out 2>&1"
