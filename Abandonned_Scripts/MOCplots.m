%% MOC PLOTTER
% Iridis runs 5km resolution
% Option 1 = Run 1
% Option 2 = Slopetest
% Option 3 = Slope33
% Option 4 = Slope7
% Option 5 = Fulltopotest
function [TS]=MOCplots(start,stop,Plot,option)
fprintf('\n dx set to 5000m')
close all
%clear all
global dx
global Z;  
global Y;
global lm;
global V;
t=10;
dx=5000;
if option ==1
%% RUN1
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Run1/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Run1/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
        fprintf(['\n',fname,' does not exist. Checking if a 5 year run was done...'])
        fname=[num2str(i),'-',num2str(i+5),'all.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning, Run1'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Run1
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Run1/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Run1
%cd /noc/users/hb1g13/Figures/Iridis4/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==2
%% Slopetest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33test/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33test/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
     if exist(fname,'file')==0
        fname=[num2str(i),'-',num2str(i+5),'all.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    [psi2]= MOC(0);%This function is quicker and matches my result
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Slope'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope33test
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope33test/glued_state_files
end
   fprintf( [num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Slope33test
%cd /noc/users/hb1g13/Figures/Iridis4/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==3
%% Slope33
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
     if exist(fname,'file')==0
        fname=[num2str(i),'-',num2str(i+5),'all.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    [psi2]= MOC(0);%This function is quicker and matches my result
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Slope with partial cells limitted to 50m min'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope33
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Slope33
%cd /noc/users/hb1g13/Figures/Iridis4/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==4
%%Slope7
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
     if exist(fname,'file')==0
        fname=[num2str(i),'-',num2str(i+5),'all.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    [psi2]=MOC(0);%This function is quicker and matches my result
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Slope with hfcamin & Adv scheme 7'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope7
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in 10 year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Slope7
%cd /noc/users/hb1g13/Figures/Iridis4/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==5
%% Fulltopotest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Fulltopotest/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Fulltopotest/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
     if exist(fname,'file')==0
        fname=[num2str(i),'-',num2str(i+5),'all.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    [psi2]= MOC(0);%This function is quicker and matches my result
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Iridis configuration, adv 33, partial cells limited'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/Fulltopotest
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Fulltopotest/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function Iridis configuration, adv 33, partial cells limited','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Fulltopotest
print('-dpng','Eulerian_Mean_timeseries_offline')
end
end