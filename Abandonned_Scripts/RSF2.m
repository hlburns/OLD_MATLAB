%% Working out the Residual Stream function from integrated meridional 
%% transport using potential temp layers
function [PsiTemp]=RSF2(Fname)
%% This funtion cacluates the redisual streamfunction in density co-orindates
% To do this simply integrate zonally and cumulatively sum up the water
% column.
% INPUTS:
% fname: A string containing the file name
%% Load Variables
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/Psi_dens/
VT=ncread(Fname,'LaVH1TH');%Layer integrated meridional transport
Y=ncread(Fname,'Yp1'); % Latitudinal distance
VT(VT==0)=NaN; % NaN for land 
%% Zonally integrate
VTdx=squeeze(nansum(VT*6666.666));
VTdx(VTdx==0)=NaN;
%% Cumulatively sum up the water column (i.e. low T to high T)
%VTdz=squeeze(nancumsum(VTdx,2)); % Looks like layers package always start
%at bottom so must sum the other way.
for i=48:-1:1
VTdz(:,i)=squeeze(nansum(VTdx(:,1:i),2));
end
%% Time average
PsiTemp=VTdz/10^6;
PsiTemptav=nanmean(VTdz,3)/10^6;
%% Plotting
Temp=csvread('Temp2');
Temp=Temp(1:48);
contourf(Y(50:301)/1000,(1:48),(PsiTemptav(50:301,1:48))',15)
set(gca,'YTicklabel',Temp(1:5:45))
title('Residual Stream function on temperature layers')
ylabel('Potential Temp ^oC')
xlabel('Meridional Distance (km)')
h=colorbar;
ylabel(h,'Transport (Sv)')
cd ~/Figures/Nautilus/Run1/
cd ~/Figures/Figs4meeting
print('Psi_res(y,t)fromT','-dpng')
%% Apply Equation of state
RhoNil=1000;
tAlpha=2*10^-4;
tref=20.0;
Rho=RhoNil*(1-(tAlpha*(Temp-tref)));
figure
[C,h]=contour(Y(50:301)/1000,(1:48),PsiTemptav(50:301,1:48)',15);
%set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2) % Labels contours
%clabel(C,h,'FontSize',6,'Color','k','Rotation',0); %Labels contour
title('Residual Stream function on density layers from temp layers')
ylabel('Potential Density kg/m^3')
xlabel('Meridional Distance (km)')
%set(gca,'YDir','reverse')
set(gca,'YTicklabel',Rho(45:-5:5))
h=colorbar;
ylabel(h,'Transport (Sv)')
print('Psi_res(y,rho)fromT','-dpng')
end
