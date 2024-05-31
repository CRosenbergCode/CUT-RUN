#!/bin/bash
# ran on CCTSI server

rmats.py --b1 b1.txt \
--b2 b2.txt \
--gtf Aedes.gtf \
--od rmatsOUT/ -t paired \
--nthread 8 \
--readLength 150 \
--tmp tmp/