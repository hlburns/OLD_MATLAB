%% Residual Stream Function in Density space
% Calls RSF.m to calculate the residual stream function over the years
% Resolution
% Option
function [TS,psi2]=Resplottermsm(start,stop,resolution,option)
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
t=5;
elseif resolution==5
    dx=5000;
    t=10;
end
fprintf(['dx set to ',num2str(dx)])
if resolution==6
    switch(option)
         case option==1
              OP='TEST';
    end
cd(['/noc/msm/scratch/students/hb1g13/Nautilus/',OP,'/Psi_dens/'])
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
    cd(['~/Figures/Nautilus/',OP,])
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd(['/noc/msm/scratch/students/hb1g13/Nautilus/',OP,'/Psi_dens/'])
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd(['~/Figures/Nautilus/',OP,])
print('-dpng','RSF_timeseries')
elseif resolution==5
      
        if option==2
                  OP='Sponge1';
        elseif option==3
                 OP='Sponge2';
        elseif option==4
                 OP='TEST';
        end
        fprintf(['\n option ',OP,' selected']);
cd(['/noc/msm/scratch/students/hb1g13/Iridis4/',OP,'/Psi_dens/'])
TS=zeros((stop-start)/t,1); 
%% Read .nc files and put into MOCcalc
t=10;
for i=start:t:stop %t year dumps  
    yr=(i-start)/t+1; %incremental steps from 1 for filling in matix
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    if exist(fname,'file')==0
       yr2=i+10;
       fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
       if exist(fname,'file')==0
          yr2=i+20;
          fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
       end
       if exist(fname,'file')==0
           yr2=i+20;
       fname=[num2str(i+10),'-',num2str(yr2),'Psi.nc'];
       end
       if exist(fname,'file')==0
           yr2=i+30;
       fname=[num2str(i+10),'-',num2str(yr2),'Psi.nc'];
       end
    end
    VT=ncread(fname,'LaVH1RHO');%Layer integrated meridional transport
    Y=ncread(fname,'Yp1'); % Latitudinal distance
    Rho=csvread('Dens');
    psi2=RSF(1);
    title('Residual Stream function on density layers, changed forcing','fontsize',12)
    cd(['~/Figures/Iridis4/',OP,])
    print(['Psires_',num2str(i),'-',num2str(yr2)],'-dpng')
    close all
    cd(['/noc/msm/scratch/students/hb1g13/Iridis4/',OP,'/Psi_dens/'])
    fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
    TS(yr)=max(max(abs(squeeze(mean(psi2,3)))));
end
figure
x=start:t:stop; %year axis in t year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd(['~/Figures/Iridis4/',OP,])
print('-dpng','RSF_timeseries')
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end