
function [TS]=MOCplottersata(start,stop,Plot,resolution,option)
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
OP='alltopo7';
%% Eularian Stream function timeseries
cd(['/noc/altix/scratch/hb1g13/Nautilus/',OP,'/grid'])
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd(['/noc/altix/scratch/hb1g13/Nautilus/',OP,'/glued_state_files'])
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
    cd(['~/Figures/Nautilus/',OP])
    print([num2str(i),'-',num2str(yr2),' MOC'],'-dpng')
    close all
    cd(['/noc/altix/scratch/hb1g13/Nautilus/',OP,'/glued_state_files'])
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
cd(['/noc/users/hb1g13/Figures/Nautilus/',OP])
print('-dpng','Eulerian_Mean_Timeseries')
elseif option==2 || option==3 && resolution==5
%% Iridis4 Alltopo7 Full diff
switch(option)
    case option==2
OP='alltopo7';
    case option==3
OP='Slope7';
end
%% Eularian Stream function timeseries
cd(['/noc/altix/scratch/hb1g13/Iridis4/',OP,'/grid'])
lm=ncread('grid.nc','HFacS');
Z=ncread('grid.nc','Zl');
Y=ncread('grid.nc','Yp1');
cd(['/noc/altix/scratch/hb1g13/Iridis4/',OP,'/glued_state_files'])
%Define outputmatrix size (no. years,t yearly avs per file) 
TS=zeros(((stop-start)/t)+5,1); 
tstep=zeros(((stop-start)/t)+5,1); 
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
    cd(['~/Figures/Iridis4/',OP,])
    print([num2str(i),'-',num2str(yr2),'MOC'],'-dpng')
    close all
    cd(['/noc/altix/scratch/hb1g13/Iridis4/',OP,'/glued_state_files'])
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
cd(['/noc/users/hb1g13/Figures/Iridis4/',OP])
print('-dpng','Eulerian_Mean_timeseries_offline')

else
    fprintf('\n Invalid combination: \n please combine option 1 with res 6km\n option 2 and 3 are set with 5km res');
end
end