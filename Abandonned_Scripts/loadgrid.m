function [GRID] = loadgrid(varargin)
%loadgrid()
%loadgrid(DIRECTORY)
%
%Reads MITgcm output files and input "data" file to create a GRID structure
%If DIRECTORY is not specified the current working directory is used.
%
%The following files are expected to be in the directory:
% data XC.* YC.* XG.* YG.* DXG.* DYG.* RAC.* hFacC.* hFacW.* hFacS.*
%
%e.g.
%>> GRID=loadgrid
%GRID =
%     drf: [50 70 100 140 190 240 290 340 390 440 490 540 590 640 690]
%     drc: [25 60 85 120 165 215 265 315 365 415 465 515 565 615 665]
%      rf: [1x16 double]
%      rc: [1x15 double]
%      xc: [90x40 double]
%      yc: [90x40 double]
%      xg: [90x40 double]
%      yg: [90x40 double]
%     dxg: [90x40 double]
%     dyg: [90x40 double]
%     rac: [90x40 double]
%   hfacc: [90x40x15 double]
%   hfacw: [90x40x15 double]
%   hfacs: [90x40x15 double]
%    mskc: [90x40x15 double]
%    mAtl: [90x40 double]
%    matl: [90x40 double]
%    mPac: [90x40 double]
%    mpac: [90x40 double]
%    msoc: [90x40 double]
%>> GRID2=loadgrid('/scratch/john/run2/');
%
%Most elements of the structure are named to corresponding to MITgcm
%variables.
%
%Lower case regional masks (matl, mpac and msoc) denote Atlantic and Pacific
%Southern Ocean regions only and are exclusive (i.e. matl does not extend
%into the Southern Ocean).
%Upper case regional masks (mAtl and mPac) denotes sectors and do include the
%Southern Ocean.
%  sum(sum( GRID.mskc(:,:,1) )) = sum( GRID.mAtl(:) + GRID.mPac(:) )
%           = sum( GRID.matl(:) + GRID.mpac(:) + GRID.msoc(:) )
%
%Written by adcroft@mit.edu, 2001
%$Header:

if nargin==0
 Dir='./';
elseif nargin==1
 Dir=[varargin{1} '/'];
else
 error('I don''t know what to do with the second argument');
end

% Extract drF from "data" file
datafile=[Dir 'data'];
fid=fopen(datafile,'r');
if fid==-1
 error(['Could not open file:' datafile ' for reading']);
end
fclose(fid);
drf=evalc([ ...
      '! grep -v ''#'' ' datafile ...
      '| awk ''/[dD[eE][lL][rRzZpP]/,/XXX/ {printf "%s",$0}'' ' ...
      '| sed ''s/[dD][eE][lL][zZrRpP][ ]*=\([0-9Ee,\. +\-]*\).*/\1/''  ' ...
      '| sed ''s/[ ]//g''  ' ...
      '| sed ''s/,/ /g''  ' ...
      ';'
     ]);
eval(['drf=[' drf '];'])
drf=drf;
drc=(drf([1 1:end-1])+drf)/2; drc(1)=drc(1)/2;

%% % Extract drC from output file
%% outputfile='output.txt';
%% drc=evalc(['!head -1000 ' outputfile ...
%%            ' | awk ''/drC/,/;/ {print $3}'' -' ...
%%            ' | egrep "e|E"' ...
%%            ' | sed ''s/,//'' ']);
%% eval(['drc=[' drc '];'])
%% drc=drc';
%% % Extract drF from output file
%% drf=evalc('!head -1000 output.txt | awk "/drF/,/;/ {print \$3}" - | egrep "e|E" | sed s/,// ');
%% eval(['drf=[' drf '];'])
%% drf=drf';

rf=-cumsum([0 drf]);
rc=-cumsum([drc]);

GRID.drf=drf;
GRID.drc=drc;
GRID.rf=rf;
GRID.rc=rc;

msg1='Error: Grid data written by the model is needed. The following files are needed: XC.* YC.* XG.* YG.* DXG.* DYG.* RAC.* hFacC.* hFacW.* hFacS.*';
msg2='Error: The appropriate model output files are not present. If the following files were not written by the model check that you have version 1.11 or greater of ini_grid.F: DXG.* DYG.* RAC.*';

xc=locrdmds([Dir 'XC'],msg1); GRID.xc=xc;
yc=locrdmds([Dir 'YC'],msg1); GRID.yc=yc;
xg=locrdmds([Dir 'XG'],msg1); GRID.xg=xg;
yg=locrdmds([Dir 'YG'],msg1); GRID.yg=yg;
dxg=locrdmds([Dir 'DXG'],msg2); GRID.dxg=dxg;
dyg=locrdmds([Dir 'DYG'],msg2); GRID.dyg=dyg;
rac=locrdmds([Dir 'RAC'],msg2); GRID.rac=rac;
hfacc=locrdmds([Dir 'hFacC'],msg1); GRID.hfacc=hfacc;
hfacw=locrdmds([Dir 'hFacW'],msg1); GRID.hfacw=hfacw;
hfacs=locrdmds([Dir 'hFacS'],msg1); GRID.hfacs=hfacs;

mskc=hfacc; mskc(find(hfacc~=0))=1; GRID.mskc=mskc;

mskc=mskc(:,:,1);
j=[];
%j=[j find( yc>-32.5 & xc>290 )'];
%j=[j find( yc>-32.5 & xc<25 )'];
j=[j find( xc>290 )'];
j=[j find( xc<25 & xc>290-360 )'];
j=[j find( yc>9 & yc<60 & (yc-9)+(xc-276)>0 )'];
j=[j find( yc>9 & yc<60 & (yc-9)+(xc-276+360)>0 & xc<30 )'];
j=[j find( yc>17 & yc<60 & xc>261 )'];
j=[j find( yc>50 & (yc-70)-(xc-270)<0 )'];
j=[j find( yc>31 & xc<38 & xc>-90 )'];
j=[j find( yc>31 & xc>360-90 )'];
j=[j find( yc>64 )'];
matl=0*mskc;
matl(j)=1;
matl=matl.*mskc;
GRID.mAtl=matl;

j=[];
j=[j find( yc>-32.5 )'];
matl=0*mskc;
matl(j)=1;
matl=GRID.mAtl.*matl.*mskc;
GRID.matl=matl;

mpac=(1-matl).*mskc;
j=[];
j=[j find( yc<65 )'];
mpac=0*mskc;
mpac(j)=1;
mpac=mpac.*mskc.*(1-GRID.mAtl);
GRID.mPac=mpac;

j=[];
j=[j find( yc>-32.5 )'];
mpac=0*mskc;
mpac(j)=1;
mpac=GRID.mPac.*mpac.*mskc;
GRID.mpac=mpac;

msoc=(1-matl).*(1-mpac).*mskc;
j=[];
j=[j find( yc<0 )'];
GRID.msoc=msoc.*mskc.*(1-mpac).*(1-matl);

function [A] = locrdmds(fname,errmsg)
[fid,msg]=fopen([fname '.001.001.meta'],'r');
if fid == -1
 A=[];
 disp(errmsg);
 error(['The error occured while trying to open ' fname '.001.001.meta'])
else
 fclose(fid);
 A=rdmds(fname);
end
