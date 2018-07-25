%% MOC offline
%This function calculates the eulerian mean overturning streamfunction
%This should be resolution, domain size and time indepedent through use of
%global variables
%INPUTS - Justs tell it whether you want it to plot of not
%OUTPUT - Timeaveraged Psi in m^3/s 
function [Psi]=MOC(Plot)
global Z
global Y
global V
global lm
global dx
%Depth integrate
dz=Z(1:end-1)-Z(2:end);
dz=[0-Z(1);dz];
lm(lm<1)=NaN; 
V=squeeze(nanmean(V(:,:,:,1:end-1),4)).*lm;
Vf=V(:,:,end:-1:1);
dzf=dz(end:-1:1);
%Zonally integrate
Vfdx=squeeze(nansum(Vf*dx));
%Depth integrate
Vfdx(Vfdx==0)=NaN;
Vdz=zeros(length(Y),length(Z)+1);
for i=1:length(Y)
    Vdz(i,2:length(Z)+1)=squeeze(Vfdx(i,:)).*dzf';% leave bottom value as zeros
end
Psi=nancumsum(Vdz,2); %Sum up the water column
Psi=Psi(:,end:-1:1); %Put back into right order 
if Plot==1   
    pcolor(Y/1000,[Z;Z(end)-250],Psi(:,1:length(Z)+1)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((Psi(:,1:length(Z)+1))));
    cmin=min(min((Psi(:,1:length(Z)+1))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title(' Eulerian mean overturning','fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
end
end