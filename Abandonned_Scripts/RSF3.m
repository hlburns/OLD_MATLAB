%% Working out RSF From Density layers
%% Working out the Residual Stream function from integrated meridional 
%% transport using potential temp layers
function [PsiRho]=RSF3(Fname,Plot)
%% This funtion cacluates the redisual streamfunction in density co-orindates
% To do this simply integrate zonally and cumulatively sum up the water
% column.
% INPUTS:
% fname: A string containing the file name
% Plot: 1 = Plot it!
%% Load Variables
%cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/Psi_dens/Rho_layers
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_flat/Psi_dens/
VT=ncread(Fname,'LaVH1RHO');%Layer integrated meridional transport
Y=ncread(Fname,'Yp1'); % Latitudinal distance
VT(VT==0)=NaN; % NaN for land 
%% Zonally integrate
VTdx=squeeze(nansum(VT*6666.666));
VTdx(VTdx==0)=NaN;
%% Cumulatively sum up the water column (i.e. low T to high T)
%VTf=flipdim(VTdx,2);
VTdz=squeeze(zeros(size(VT(1,:,:,:))));
for i=length(VTdx(1,:,1)):-1:1
VTdz(:,i,:)=squeeze(nansum(VTdx(:,1:i,:),2));
end
%VTdz=flipdim(VTdz,2);
VTdz(VTdz==0)=NaN;
%% Time average
PsiRho=VTdz/10^6;
PsiRhotav=nanmean(VTdz,3)/10^6;
%% Plotting
if Plot==1
Rho=csvread('Dens50');
Rho=Rho(1:58);
contourf(Y(1:301)/1000,(1:58),(PsiRhotav(1:301,1:58))',10)
shading flat
cmax=max(max((PsiRhotav)));
cmin=min(min((PsiRhotav)));
colormap(b2r(cmin,cmax))
set(gca,'YTicklabel',Rho(1:5:58))
title('Residual Stream function on density layers')
ylabel('Potential Density kg/m^3')
set(gca,'YDir','reverse')
xlabel('Meridional Distance (km)')
h=colorbar;
ylabel(h,'Transport (Sv)')
%cd ~/Figures/Nautilus/Run1
cd ~/Figures/Nautilus/Flat/
%cd ~/Figures/Figs4meeting
print('Psi_res(y,rho)dlayersflat2','-dpng')
end
end
