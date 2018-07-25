%% Residual Stream Function in Density space
% Calls RSF.m to calculate the residual stream function over the years
% Resolution
% Option
function [TS,psi2]=Resplot_it(start,stop,t,resolution,option)
close all
if nargin ~= 5 ;
   disp('input error');
   fprintf('input 5 variables: \n Start year, Stop year,\n t=timesteps,\n Resolution 5 or 6km,\n Options:\n 1= Nautilus_fulltopo_nodiff, \n 2 = Iridis fulltopo full diff,\n Iridis Slope7 ')
end
global dx  
global Y;
global VT;
global Rho;
if resolution==6
dx=6666.6666667;
end
fprintf(['dx set to',num2str(dx)])
if option==1 && resolution==6
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/Psi_dens/
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(0);
    Psi(:,:,yr)=psi2;
    fprintf(['\n',num2str(yr),' out of ',num2str((start-stop)/t)])
    yr=yr+1;
end
    psi2=squeeze(nanmean(Psi,3));
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    title('Residual Stream function on density layers, Alltopo7 nodiff','fontsize',12)
    cd ~/Figures/Nautilus/alltopo7
    print(['Psires_',num2str(Start),'-',num2str(stop)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/Psi_dens/
 
elseif option==2 && resolution==5
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/Psi_dens/
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(0);
    Psi(:,:,yr)=psi2;
    fprintf(['\n',num2str(yr),' out of ',num2str((start-stop)/t)])
    yr=yr+1;
end
    psi2=squeeze(nanmean(Psi,3));
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    title('Residual Stream function on density layers,\n Alltopo7 fulldiff','fontsize',12)
    cd ~/Figures/Iridis4/alltopo7/
    print(['Psires',num2str(start),'-',num2str(stop)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/Psi_dens/
elseif option==3 && resolution==5
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(0);
    Psi(:,:,yr)=psi2;
    fprintf(['\n',num2str(yr),' out of ',num2str((start-stop)/t)])
    yr=yr+1;
end
    psi2=squeeze(nanmean(Psi,3));
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    
    title('Residual Stream function on density layers, Slope7 fulldiff','fontsize',12)
    cd ~/Figures/Iridis4/Slope7/
    print(['Psires',num2str(start),'-',num2str(stop)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end