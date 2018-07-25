function [TS]=ThetaTS(start,stop)
%% Eularian Stream function timeseries
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo/grid
global lm;
lm=ncread('grid.nc','HFacC');
cd /noc/altix/scratch/hb1g13/Nautilus/alltopo/glued_state_files
%Define outputmatrix size (no. years,5 yearly avs per file) 
TS=zeros((stop-start)/5,1); 
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'all.nc'];
    T=ncread(fname,'Temp');
    Temp=squeeze(mean(T,4)).*lm;
    Theta=squeeze(mean(squeeze(mean(squeeze(mean(Temp))))));
   TS(yr)=Theta;
end
figure
x=start:5:stop; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Mean Basin Temperature','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Theta ^oC','fontsize',12)
%cd /noc/users/hb1g13/Figures/Nautilus/Flat
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Mean_Theta_timeseries')
end