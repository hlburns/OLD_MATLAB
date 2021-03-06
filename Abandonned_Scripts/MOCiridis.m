function [TS]=EulTSoffline2(start,stop,Plot,option)
if option ==1
%% Run1
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 6666.666;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    %V=ncread(fname,'V');
    %V=V(:,:,end:-1:1,:);
    %vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    %[psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
    [~,psi]=MOCoffline(fname,0);
if Plot==1
    %psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    psi2=psi;
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((psi2(:,1:25))));
    cmin=min(min((psi2(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/Adv/
    %cd ~/Figures/Nautilus/Flat
    print([num2str(i),'-',num2str(yr2),'SFoffiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Adv
print('-dpng','Eulerian_Mean_timeseries_offline')
end
if option==2
%% Slopetest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 6666.666;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    %V=ncread(fname,'V');
    %V=V(:,:,end:-1:1,:);
    %vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    %[psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
    [~,psi]=MOCoffline(fname,0);
if Plot==1
    %psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    psi2=psi;
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((psi2(:,1:25))));
    cmin=min(min((psi2(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/Adv/
    %cd ~/Figures/Nautilus/Flat
    print([num2str(i),'-',num2str(yr2),'SFoffiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Adv
print('-dpng','Eulerian_Mean_timeseries_offline')
end
if option==3
%% Slope33
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 6666.666;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    %V=ncread(fname,'V');
    %V=V(:,:,end:-1:1,:);
    %vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    %[psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
    [~,psi]=MOCoffline(fname,0);
if Plot==1
    %psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    psi2=psi;
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((psi2(:,1:25))));
    cmin=min(min((psi2(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/Adv/
    %cd ~/Figures/Nautilus/Flat
    print([num2str(i),'-',num2str(yr2),'SFoffiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Adv
print('-dpng','Eulerian_Mean_timeseries_offline')
end
if option==3
%% Slope7
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 6666.666;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    %V=ncread(fname,'V');
    %V=V(:,:,end:-1:1,:);
    %vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    %[psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
    [~,psi]=MOCoffline(fname,0);
if Plot==1
    %psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    psi2=psi;
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((psi2(:,1:25))));
    cmin=min(min((psi2(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/Adv/
    %cd ~/Figures/Nautilus/Flat
    print([num2str(i),'-',num2str(yr2),'SFoffiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Adv
print('-dpng','Eulerian_Mean_timeseries_offline')
end
if option==4
%% Alltopotest
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/grid
global Z; 
global Zl; 
global dz;
global dx;
global Y;
global lm;
lm=ncread('grid.nc','HFacS');
dx = 6666.666;
Zl=ncread('grid.nc','Zl');
Z=ncread('grid.nc','Z');
Y=ncread('grid.nc','Yp1');
dz=Zl(1:23)-Zl(2:24);
dz=[0-Z(1);dz;200];
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    %V=ncread(fname,'V');
    %V=V(:,:,end:-1:1,:);
    %vmask=change(lm,'>',0,1); %This is now EXACLY what the script calls for
    %vmask=vmask(:,:,end:-1:1);
    %vmask(~isnan(vmask))=1;
    %[psi,~]= mit_overturning(V,vmask,dx,dz,1);%This function is quicker and matches my result
    [~,psi]=MOCoffline(fname,0);
if Plot==1
    %psi2=squeeze(nanmean(psi(:,:,1:60),3)); %bin the last time step
    psi2=psi;
    contourf(Y/1000,[Zl;-4000],psi2(:,1:25)');  %can add ,15 to add more contours
    shading flat
    cmax=max(max((psi2(:,1:25))));
    cmin=min(min((psi2(:,1:25))));
    colormap(b2r(cmin,cmax))
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/Adv/
    %cd ~/Figures/Nautilus/Flat
    print([num2str(i),'-',num2str(yr2),'SFoffiline'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel_adv/glued_state_files
end
   TS(yr)=max(max(abs(squeeze(nanmean(psi,3)))));
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Adv
print('-dpng','Eulerian_Mean_timeseries_offline')
end
end
