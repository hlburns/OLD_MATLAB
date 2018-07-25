%% MOC 20 year av
% This function will plot a set of MOCs and produce a timeseries.
% This function calls MOC.m
% Please set resolution
% Option 1 = Nautilus alltopo_nodiff
% Option 2 = Iridis alltopo_fulldiff
% Option 3 = Iridis Slope 7 (full diff)
function [Psi]=MOCplot_it(start,stop,t,resolution,option)
close all
if nargin ~= 5 ;
   disp('input error');
   fprintf('input 4 variables: \n Start year, Stop year,\n t=timestep in years,\n Resolution 5 or 6km,\n Options:\n 1= Nautilus_fulltopo_nodiff, \n 2 = Iridis fulltopo full diff,\n Iridis Slope7 ')
end
global dx
global Z;  
global Y;
global lm;
global V;
dx=5000;
if resolution==6
dx=6666.6666667;
end
fprintf(['dx set to ',num2str(dx)])
if option ==1 && resolution==6
%% Nautilus alltopo_nodiff
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/glued_state_files
%% Read .nc files and put into MOCcalc
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    [psi2]= MOC(0);
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
    title([num2str(start),'-',num2str(stop),' Eulerian mean overturning, Scheme 7,tref=0, Hfacmin set, nodiffusion'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/alltopo7
    print([num2str(start),'-',num2str(stop),' MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/glued_state_files

elseif option==2 && resolution==5
%% Iridis4 Alltopo7 Full diff
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
%% Read .nc files and put into MOCcalc
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    [psi2]= MOC(0);
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
    title([num2str(start),'-',num2str(stop),' Eulerian mean overturning, Fulltopo, HFaCmin set, Adv scheme 7, full diffusion '],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/alltopo7
    print([num2str(start),'-',num2str(stop),'MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/glued_state_files

elseif option==3 && resolution==5
%% Slope7
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
%% Read .nc files and put into MOCcalc
yr=1;
Psi=zeros(length(Y),25,(start-stop)/t);
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    [psi2]= MOC(0);
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
    title([num2str(start),'-',num2str(stop),' Eulerian mean overturning Slope only, Full diffusion'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/Slope7
    print([num2str(start),'-',num2str(stop),'MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end