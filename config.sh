#!/bin/bash
ssh $1 'sudo apt-get update'
ssh $1 'sudo apt-get -y install python3-pip'
scp dependencies.txt $1:~
ssh $1 'sudo pip3 install -r dependencies.txt'
ssh $1 'sudo pip3 install cvxpy'
ssh $1 'sudo apt-get -y install python3-tk'
ssh $1 'mkdir sat_figures saturated_data recon_sat_figures network_recon_sat_figures'
echo 'Complete'

