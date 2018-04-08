#!/bin/bash
scp -r $1:~/saturated_data/ saturated_data
scp -r $1:~/sat_figures/ sat_figures/
scp -r $1:~/recon_sat_figures/ recon_sat_figures/
scp -r $1:~/network_recon_sat_figures/ network_recon_sat_figures/
