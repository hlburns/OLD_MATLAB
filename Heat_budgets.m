%% Heat Budget Stuff
% The data:
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/Heat_Bugets
fname='Heat.nc';
fname2='Heat_Surface.nc';
Tend=ncread(fname,'TOTTTEND'); %X,Y,Z (degc/day)
Avx1=ncread(fname,'ADVx_TH'); %Xp1,Y,Z (degC.m^3/s)
Avy1=ncread(fname,'ADVy_TH'); %X,Yp1,Z (degC.m^3/s)
AVz=ncread(fname,'ADVr_TH'); %X,Y,Z  (degC.m^3/s)
Dfx1=ncread(fname,'DFxE_TH'); %Xp1,Y,Z (degC.m^3/s)
Dfy1=ncread(fname,'DFyE_TH'); %X,Yp1,Z (degC.m^3/s)
Dfz=ncread(fname,'DFrE_TH'); %X,Y,Z (degC.m^3/s)
TFLUX=ncread(fname2,'TFLUX'); %X,Y (W/m^2)
% Normally things are in W so convert to Watts
% Grid file lm
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/grid
lm= ncread('grid.nc','HFacC'); % % box open at center of box
%% Regrid all to same points all X,Y,Z
% NB Xp1 is West face of cell
%    Yp1 is South face of cell
%% Create dV boxes and dA boxes for surface
%% Convert to Watts

