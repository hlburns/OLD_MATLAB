%% MOC PLOTTER
% This function will plot a set of MOCs and produce a timeseries.
% This function calls MOC.m
% Please set resolution
% Option 1 = Nautilus TEST
% Option 2 = Iridis Sponge1
% Option 3 = Iridis Sponge2
% Option 4 = Iridis TEST
function [TS]=MOCplottermsm(start,stop,resolution,option)
close all
if nargin ~= 4 ;
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
t=5;
elseif resolution==5
    dx=5000;
    t=10;
end
fprintf(['dx set to ',num2str(dx)])
if resolution==6
%% Changed Forcing
switch(option)
    case option==1
         OP='TEST';
end
%% Eularian Stream function timeseries
cd(['/noc/msm/scratch/students/hb1g13/Nautilus/',OP,'/grid'])
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd(['/noc/msm/scratch/students/hb1g13/Nautilus/',OP,'/glued_state_files'])
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
    cd(['~/Figures/Nautilus/',OP])
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd(['/noc/msm/scratch/students/hb1g13/Nautilus/',OP,'/glued_state_files'])
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
cd(['/noc/users/hb1g13/Figures/Nautilus/',OP])
print('-dpng','Eulerian_Mean_timeseries')
elseif resolution==5
%% Sponge1 Sponge2 TEST
if option==2
         OP='Sponge1';
elseif option==3
         OP='Sponge2';
elseif option==4
         OP='TEST';
end
%% Eularian Stream function timeseries
cd(['/noc/msm/scratch/students/hb1g13/Iridis4/',OP,'/grid'])
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd(['/noc/msm/scratch/students/hb1g13/Iridis4/',OP,'/glued_state_files'])
%Define outputmatrix size (no. years,10 yearly avs per file) 
TS=zeros(((stop-start)/t)+5,1); 
tstep=zeros(((stop-start)/t)+5,1);
%% Read .nc files and put into MOCcalc
yr=1;
for i=start:t:stop %t year dumps
    t=10;
    yr2=i+t; %t year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    if exist(fname,'file')==0
       yr2=i+5;
       fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       if exist(fname,'file')==0
          yr2=i+20;
          fname=[num2str(i),'-',num2str(yr2),'all.nc'];
       end
       if exist(fname,'file')==0
           yr2=i+20;
       fname=[num2str(i+10),'-',num2str(yr2),'all.nc'];
       end
    end
    V=ncread(fname,'V');
    if length(V(1,1,1,:)) <120
    t=5;
    elseif length(V(1,1,1,:))>200
        t=20;
    else
      t=10;
    end
    yr=yr+1;
    [psi2]= MOC(0);
    contourf(Y/1000,[Z;Z(end)-250],psi2'/10^6,15);  %can add ,15 to add more contours
    cmax=max(max((psi2(:,1:end-1)/10^6)));
    cmin=min(min((psi2(:,1:end-1)/10^6)));
    colormap(b2r(cmin,cmax))
    shading flat
    xlabel('Meridional distance (km)','fontsize',12)
    ylabel('Depth (m)','fontsize',12)
    title([num2str(i),'-',num2str(yr2),' Eulerian mean overturning ',OP,''],'fontsize',12)
    h=colorbar;
    ylabel(h,'Transport (Sv)','fontsize',12)
    cd(['~/Figures/Iridis4/',OP])
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd(['/noc/msm/scratch/students/hb1g13/Iridis4/',OP,'/glued_state_files'])
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
cd(['/noc/users/hb1g13/Figures/Iridis4/',OP])
print('-dpng','Eulerian_Mean_timeseries')
else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end