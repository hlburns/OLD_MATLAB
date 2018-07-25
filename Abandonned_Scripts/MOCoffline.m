%% MOC offline
function [maxosf, Psi]=MOCoffline(fname,Plot)
global Z; 
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files/
V=ncread(fname,'V');
Y=ncread(fname,'Yp1');
%X=ncread('fname.nc','X');
Z=ncread(fname,'Zl');
%Depth integrate
dz=Z(1:23)-Z(2:24);
dz=[0-Z(1);dz];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid/
lm=ncread('grid.nc','HFacS');
lm(lm<1)=NaN; 
V=squeeze(nanmean(V,4)).*lm;
%Vf=flipdim(V,3);%Flip to integrate bottom up
%dzf=flipdim(dz,1);
Vf=V(:,:,end:-1:1);
dzf=dz(end:-1:1);
%Zonally integrate
dx=6666.666; 
Vfdx=squeeze(nansum(Vf*dx));
%Depth integrate
Vfdx(Vfdx==0)=NaN;
Vdz=zeros(301,25);
for i=1:301
    Vdz(i,2:25)=squeeze(Vfdx(i,:)).*dzf';
end
Psi=nancumsum(Vdz,2)/10^6;
Psi=Psi(:,end:-1:1);
maxosf=max(max(abs(Psi)));
if Plot==1   
    contourf(Y/1000,[Z;-4000],Psi(:,1:25)');  %can add ,15 to add more contours
    cmax=max(max((Psi(:,1:25))));
    cmin=min(min((Psi(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([fname,' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
end
end