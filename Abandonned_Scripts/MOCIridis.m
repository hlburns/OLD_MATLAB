function [TS]=EulTSoffline(start,stop,Plot,option)
if option ==1
%% RUN1
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Run1/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 5000;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/Iridis4/Run1/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/10,1); 
%% Read .nc files and put into MOCcalc
for i=start:10:stop %10 year dumps  
    yr=(i-start)/10+1; %incremental steps from 1 for filling in matix
    yr2=i+10; %10 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    V=V(:,:,end:-1:1,:);
    vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %%vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    [psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
if Plot==1
    psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)'/10^6,15,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    %shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Run1
    print([num2str(i),'-',num2str(yr2),'SF1offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Run1/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:10:stop; %year axis in 10 year time steps and then yearly steps
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
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 5000;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33test/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/10,1); 
%% Read .nc files and put into MOCcalc
for i=start:10:stop %10 year dumps  
    yr=(i-start)/10+1; %incremental steps from 1 for filling in matix
    yr2=i+10; %10 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    V=V(:,:,end:-1:1,:);
    vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    [psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
if Plot==1
    psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    %shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope33test
    print([num2str(i),'-',num2str(yr2),'SF1offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope33test/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:10:stop; %year axis in 10 year time steps and then yearly steps
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
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 5000;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/10,1); 
%% Read .nc files and put into MOCcalc
for i=start:10:stop %10 year dumps  
    yr=(i-start)/10+1; %incremental steps from 1 for filling in matix
    yr2=i+10; %10 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    V=V(:,:,end:-1:1,:);
    vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    [psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
if Plot==1
    psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    %shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope33
    print([num2str(i),'-',num2str(yr2),'SF1offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:10:stop; %year axis in 10 year time steps and then yearly steps
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
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 5000;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/10,1); 
%% Read .nc files and put into MOCcalc
for i=start:10:stop %10 year dumps  
    yr=(i-start)/10+1; %incremental steps from 1 for filling in matix
    yr2=i+10; %10 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    V=V(:,:,end:-1:1,:);
    vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    [psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
if Plot==1
    psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    %shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    %cd ~/Figures/Iridis4/Run1/
    cd ~/Figures/Iridis4/Slope7
    print([num2str(i),'-',num2str(yr2),'SF1offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:10:stop; %year axis in 10 year time steps and then yearly steps
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
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 5000;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/Iridis4/Fulltopotest/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/10,1); 
%% Read .nc files and put into MOCcalc
for i=start:10:stop %10 year dumps  
    yr=(i-start)/10+1; %incremental steps from 1 for filling in matix
    yr2=i+10; %10 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    V=V(:,:,end:-1:1,:);
    vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    [psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
if Plot==1
    psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    %shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Iridis configuration, adv 33, artial cells limited'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/Fulltopotest
    print([num2str(i),'-',num2str(yr2),'SF1offiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Fulltopotest/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:10:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function Iridis configuration, adv 33, partial cells limited','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Fulltopotest
print('-dpng','Eulerian_Mean_timeseries_offline')
end
end
