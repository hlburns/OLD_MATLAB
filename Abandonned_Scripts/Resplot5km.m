%% Residual Stream Function in Density space
function [TS,psi2]=Resplot5km(start,stop,option)
fprintf('dx set to 5km')
close all
%clear all
global dx  
global Y;
global VT;
global Rho;
t=5;
dx=5000;
if option==1
cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/Psi_dens/
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
    cd ~/Figures/Iridis4/Slope33/
    print(['Psires_Slope33',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope33/Psi_dens/
    fprintf([num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/Slope33/
print('-dpng','RSF_timeseries_offline')
end
if option==2
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    if exist(fname,'file')==0
        fname=[num2str(i),'-',num2str(i+5),'Psi.nc'];
        if exist(fname,'file')==0
        fname=[num2str(i-5),'-',num2str(i+5),'Psi.nc'];
        end
    end
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    cd ~/Figures/Iridis4/Slope7/
    print(['Psires_Slope7,',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/Psi_dens/
    fprintf([num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd ~/Figures/Iridis4/Slope7/
print('-dpng','RSF_timeseries_offline')
end
end