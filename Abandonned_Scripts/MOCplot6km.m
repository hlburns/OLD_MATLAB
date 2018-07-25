%% MOC PLOTTER
% Nautilus runs 6.6666km resolution
% Option 1 = Run 1
% Option 2 = Alltopo test
% Option 3 = Flat bottom
% Optopn 4 = Alltopo Av Scheme 7
function [TS]=MOCplot6km(start,stop,Plot,option)
fprintf('dx set to 6666.666m')
close all
%clear all
global dx
global Z;  
global Y;
global lm;
global V;
t=5;
dx=6666.6666667;
if option<1 || option>3
    fprintf('\n This is the Nautilus 6km res plotter: Option 1 = Run1, Option 2 = alltopotest and Option 3 = Flat bottom')
end
if option ==1
%% RUN1
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/nchannel/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/nchannel/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
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
    %cd ~/Figures/Nautilus/Run1/
    cd ~/Figures/Nautilus/Run1
    print([num2str(i),'-',num2str(yr2),'offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/nchannel/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
%cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==2
%% Alltopotest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
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
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Scheme 33, tref=0, Hfacmin set'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Nautilus/Run1/
    cd ~/Figures/Nautilus/alltopo
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/alltopo
%cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==3
%% nchannel_flat
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/nchannel_flat/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/nchannel_flat/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
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
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Flat bottom'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Nautilus/Run1/
    cd ~/Figures/Nautilus/nchannel_flat
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/nchannel_flat/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in 10 year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/nchannel_flat
%cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==4
%%All topo adv scheme 7
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
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
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Alltopo  with hfcamin & Adv scheme 7'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Nautilus/Run1/
    cd ~/Figures/Nautilus/alltopo7
    print([num2str(i),'-',num2str(yr2),'Stream function'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/glued_state_files
end
   fprintf([num2str(i),'-',num2str(yr2),' done'])
   TS(yr)=max(max(abs(psi2)));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/alltopo7
%cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Eulerian_Mean_timeseries_offline1')
end
if option==5
%% Fulltopotest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/Fulltopotest/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Nautilus/Fulltopotest/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
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
    cd ~/Figures/Nautilus/Fulltopotest
    print([num2str(i),'-',num2str(yr2),'SF2offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/Fulltopotest/glued_state_files
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
cd /noc/users/hb1g13/Figures/Nautilus/Fulltopotest
print('-dpng','Eulerian_Mean_timeseries_offline')
end
end
