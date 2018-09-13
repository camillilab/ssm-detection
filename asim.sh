#!/bin/sh

python3 art_sim.py -p HS25 -n 10000000 -R ./ssm-detection/scripts/art_refs/ -s 10 -l 150 -i 250 | tee artsim.log
