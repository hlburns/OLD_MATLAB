%% Residual Stream Function in Density space
% Calls RSF.m to calculate the residual stream function over the years
% Resolution
% Option
function [TS,psi2]=Resplotter(start,stop,resolution,option)
close all
if nargin ~= 4 ;
   disp('input error');
   fprintf('input 5 variables: \n Start year, Stop year,\n Plot=1 or no Plot=0,\n Resolution 5 or 6km,\n Options:\n 1= Nautilus_fulltopo_nodiff, \n 2 = Iridis fulltopo full diff,\n Iridis Slope7 ')
end
global dx  
global Y;
global VT;
global Rho;
if resolution==6
dx=6666.6666667;
t=4;
elseif resolution==5
    dx=5000;
    t=10;
end
fprintf(['dx set to',num2str(dx)])
if option==1 && resolution==6
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/Psi_dens/
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    title('Residual Stream function on density layers, Alltopo7 nodiff','fontsize',12)
    cd ~/Figures/Nautilus/alltopo7
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Nautilus/alltopo7
print('-dpng','RSF_timeseries')
elseif option==2 && resolution==5
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/Psi_dens/
TS=zeros((stop-start)/t,1); 
tstep=zeros((stop-start)/t+5,1); 
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    t=10;
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    if exist(fname,'file')==0
       yr2=yr2-5;
       t=5;
       fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       if exist(fname,'file')==0
          yr2=yr2+10;
          fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       end
    end
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    yr=yr+1;
    psi2=RSF(1);
    title('Residual Stream function on density layers,\n Alltopo7 fulldiff','fontsize',12)
    cd ~/Figures/Iridis4/alltopo7/
    print(['Psires',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
    tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in t year time steps and then yearly steps
x(x==0)=NaN;
TS(TS==0)=NaN;
if length(TS)==length(x)
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/alltopo7/
print('-dpng','RSF_timeseries_offline')
end
elseif option==3 && resolution==5
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
TS=zeros((stop-start)/t,1); 
tstep=zeros((stop-start)/t+5,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    if exist(fname,'file')==0
       yr2=yr2-5;
       fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
       if exist(fname,'file')==0
          yr2=yr2+10;
          fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
       end
    end
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    if length(VT(1,1,1,:))>5 &&length(VT(1,1,1,:)) <10
        t=5;
    else 
        t=10;
    end
    psi2=RSF(1);
    title('Residual Stream function on density layers, Slope7 fulldiff','fontsize',12)
    cd ~/Figures/Iridis4/Slope7/
    print(['Psires',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
    tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in t year time steps and then yearly steps
x(x==0)=NaN;
TS(TS==0)=NaN;
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/Slope7/
print('-dpng','RSF_timeseries_offline')
elseif option==4 && resolution==6
cd /noc/msm/scratch/students/hb1g13/Nautilus/TEST/Psi_dens/
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
t=5;
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    title('Residual Stream function on density layers, changed forcing','fontsize',12)
    cd ~/Figures/Nautilus/TEST
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Nautilus/TEST/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Nautilus/TEST
print('-dpng','RSF_timeseries')
elseif option==5 && resolution==5
cd /noc/msm/scratch/students/hb1g13/Iridis4/TEST/Psi_dens/
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
t=10;
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    title('Residual Stream function on density layers, changed forcing','fontsize',12)
    cd ~/Figures/Iridis4/TEST
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Iridis4/TEST/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/TEST
print('-dpng','RSF_timeseries')
elseif option==6 && resolution==5
cd /noc/msm/scratch/students/hb1g13/Iridis4/Sponge1/Psi_dens/
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
t=10;
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    title('Residual Stream function on density layers, changed forcing','fontsize',12)
    cd ~/Figures/Iridis4/Sponge1
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Iridis4/Sponge1/Psi_dens/
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/Sponge1
print('-dpng','RSF_timeseries')
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end
