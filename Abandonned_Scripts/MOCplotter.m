%% MOC PLOTTER
% This function will plot a set of MOCs and produce a timeseries.
% This function calls MOC.m
% Please set resolution
% Option 1 = Nautilus alltopo_nodiff
% Option 2 = Iridis alltopo_fulldiff
% Option 3 = Iridis Slope 7 (full diff)
% Option 4 = Nautilus TEST
% Option 5 = Iridis Sponge1
% Option 6 = Iridis TEST
function [TS]=MOCplotter(start,stop,Plot,resolution,option)
close all
if nargin ~= 5 ;
   disp('input error');
   fprintf('input 5 variables: \n Start year, Stop year,\n Plot=1 or no Plot=0,\n Resolution 5 or 6km,\n Options:\n 1= Nautilus_fulltopo_nodiff, \n 2 = Iridis fulltopo full diff,\n Iridis Slope7 ')
end
global dx
global Z;  
global Y;
global lm;
global V;
if resolution==6
dx=6666.6666667;
t=4;
elseif resolution==5
    dx=5000;
    t=10;
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
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning, Scheme 7,tref=0, Hfacmin set, nodiffusion'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/alltopo7
    print([num2str(i),'-',num2str(yr2),' MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Nautilus/alltopo7/glued_state_files
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
cd /noc/users/hb1g13/Figures/Nautilus/alltopo7
print('-dpng','Eulerian_Mean_Timeseries')
elseif option==2 && resolution==5
%% Iridis4 Alltopo7 Full diff
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/glued_state_files
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros((stop-start)/t+5,1); 
tstep=zeros((stop-start)/t+5,1); 
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
        yr2=yr2-5;
        fname=[num2str(i),'-',num2str(yr2),'all.nc'];
        if exist(fname,'file')==0
            yr2=yr2+10;
            fname=[num2str(i),'-',num2str(yr2),'all.nc'];
        end
    end
    V=ncread(fname,'V');
    t=10;
    if length(V(1,1,1,:)) <121
        t=5;
    end
    yr=yr+1;
     % The time steps may change 5,10 or 20
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning, Fulltopo, HFaCmin set, Adv scheme 7, full diffusion '],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/alltopo7
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/alltopo7/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
   tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in t year time steps and then yearly steps
x(x==0)=NaN;
TS(TS==0)=NaN;
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/alltopo7
print('-dpng','Eulerian_Mean_timeseries_offline')
elseif option==3 && resolution==5
%% Slope7
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t+5,1); 
tstep=zeros((stop-start)/t+5,1);
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
       yr2=yr2-5;
       fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       if exist(fname,'file')==0
          yr2=yr2+10;
          fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       end
    end
    V=ncread(fname,'V');
    if length(V(1,1,1,:)) <120
    t=5;
    else
      t=10;
    end
    yr=yr+1;
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Slope only, Full diffusion'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/Slope7
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd /noc/altix/scratch/hb1g13/Iridis4/Slope7/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
   tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in 10 year time steps and then yearly steps
x(x==0)=NaN; % Preallocation is difficult as t varies!
TS(TS==0)=NaN;
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Slope7
print('-dpng','Eulerian_Mean_timeseries')
elseif option==4 && resolution==6
%% Changed Forcing
%% Eularian Stream function timeseries
cd /noc/msm/scratch/students/hb1g13/Nautilus/TEST/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/msm/scratch/students/hb1g13/Nautilus/TEST/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t+5,1); 
tstep=zeros((stop-start)/t+5,1);
t=5;
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    V=ncread(fname,'V');
    yr=yr+1;
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-4000],psi2(:,1:25)'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:25)/10^6)));
    cmin=min(min((psi2(:,1:25)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Slope only, Full diffusion'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Nautilus/TEST
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Nautilus/TEST/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
   tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in 10 year time steps and then yearly steps
x(x==0)=NaN; % Preallocation is difficult as t varies!
TS(TS==0)=NaN;
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/TEST
print('-dpng','Eulerian_Mean_timeseries')
elseif option==5 && resolution==5
%% Sponge1
%% Eularian Stream function timeseries
cd /noc/msm/scratch/students/hb1g13/Iridis4/Sponge1/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/msm/scratch/students/hb1g13/Iridis4/Sponge1/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t+5,1); 
tstep=zeros((stop-start)/t+5,1);
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
       yr2=yr2-5;
       fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       if exist(fname,'file')==0
          yr2=yr2+10;
          fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       end
    end
    V=ncread(fname,'V');
    if length(V(1,1,1,:)) <120
    t=5;
    else
      t=10;
    end
    yr=yr+1;
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-3000],psi2'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:30)/10^6)));
    cmin=min(min((psi2(:,1:30)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning Sponge1'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/Sponge1
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Iridis4/Sponge1/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
   tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in 10 year time steps and then yearly steps
x(x==0)=NaN; % Preallocation is difficult as t varies!
TS(TS==0)=NaN;
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/Sponge1
print('-dpng','Eulerian_Mean_timeseries')
elseif option==6 && resolution==5
%% TEST
%% Eularian Stream function timeseries
cd /noc/msm/scratch/students/hb1g13/Iridis4/TEST/grid
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd /noc/msm/scratch/students/hb1g13/Iridis4/TEST/glued_state_files
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros((stop-start)/t+5,1); 
tstep=zeros((stop-start)/t+5,1);
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps  
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
       yr2=yr2-5;
       fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       if exist(fname,'file')==0
          yr2=yr2+10;
          fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       end
    end
    V=ncread(fname,'V');
    if length(V(1,1,1,:)) <120
    t=5;
    else
      t=10;
    end
    yr=yr+1;
    [psi2]= MOC(0);
if Plot==1
    contourf(Y/1000,[Z;-3000],psi2'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:24)/10^6)));
    cmin=min(min((psi2(:,1:24)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning TEST'],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd ~/Figures/Iridis4/TEST
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd /noc/msm/scratch/students/hb1g13/Iridis4/TEST/glued_state_files
end
   fprintf(['\n',num2str(i),'-',num2str(yr2),' done']);
   TS(yr)=max(max(abs(psi2)));
   tstep(yr)=t;
end
figure
x=cumsum(tstep)+start; %year axis in 10 year time steps and then yearly steps
x(x==0)=NaN; % Preallocation is difficult as t varies!
TS(TS==0)=NaN;
plot(x,TS/10^6,'k','linewidth',1.8);
title('Time series of Eularian Stream function offline','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Iridis4/TEST
print('-dpng','Eulerian_Mean_timeseries')
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end
