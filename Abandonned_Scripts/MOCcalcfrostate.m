%% MOC calc offline
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/glued_state_files/
V=ncread('220-225all.nc','V');
Y=ncread('220-225all.nc','Yp1');
X=ncread('220-225all.nc','X');
Z=ncread('220-225all.nc','Z');
%Depth integrate
dz=Z(1:23)-Z(2:24);
dz=[0-Z(1);dz];
V(V==0)=NaN;
Vf=flipdim(V,3);
dzf=flipdim(dz,1);
%timeaverage
Vf=squeeze(nanmean(Vf,4));
%Zonally integrated
dx=6666.666; %m
Vfdx=squeeze(nansum(Vf*dx));
%Depth integrate
Vfdx(Vfdx==0)=NaN;
Vdz=zeros(301,24);
for i=1:301
    Vdz(i,:)=squeeze(Vfdx(i,:)).*dzf';
end
Psi=nancumsum(Vdz,2)/10^6;
Psi=flipdim(Psi,2);
maxosf=max(max(abs(Psi)));