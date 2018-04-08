#!/bin/bash
ssh $1 'sudo apt-get update'
ssh $1 'sudo apt-get -y install python3-pip'
scp dependencies.txt $1:~
ssh $1 'sudo pip3 install -r dependencies.txt'
ssh $1 'sudo pip3 install cvxpy'
ssh $1 'sudo apt-get -y install python3-tk'

echo 'Complete'

