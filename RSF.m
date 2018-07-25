%% Residual Stream function Plotter
function [PsiRho]=RSF(Plot)
%% This funtion cacluates the redisual streamfunction in density co-orindates
% To do this simply integrate zonally and cumulatively sum up the water
% column.
% INPUTS:
% fname: A string containing the file name
% Plot: 1 = Plot it!
global VT
global Y
global Rho
global dx
VT(VT==0)=NaN; % NaN for land 
%% Zonally integrate
VTdx=squeeze(nansum(VT*dx));
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
Rho=Rho(1:length(PsiRhotav(1,:)));
contourf(Y/1000,(1:length(Rho)),PsiRhotav',10)
shading flat
cmax=max(max((PsiRhotav)));
cmin=min(min((PsiRhotav)));
colormap(b2r(cmin,cmax))
set(gca,'YTicklabel',Rho(1:5:end))
title('Residual Stream function on density layers','fontsize',12)
ylabel('Potential Density kg/m^3','fontsize',12)
set(gca,'YDir','reverse')
xlabel('Meridional Distance (km)','fontsize',12)
h=colorbar;
ylabel(h,'Transport (Sv)','fontsize',12)
end
end