cd /fibre/hb1g13/MITgcm
tic
V=ncread('30-35all.nc','V');
toc
clear all
cd /sata/scratch/hb1g13/Iridis4/Storecupboard/Slope33/glued_state_files
tic 
V=ncread('30-35all.nc','V');
toc