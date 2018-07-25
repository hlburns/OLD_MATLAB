function [TS]=ResTS(start,stop,option)
%% Eularian Stream function timeserie
%% Read .nc files and put into MOCcalc
if option ==1
for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    [PsiTemp]=RSF2(fname); %Call my MOCcalc function
    Psires1(1:5)=max(max(abs(PsiTemp(50:301,:,:))));
    Psires2(1:5)=max(max(abs(PsiTemp(150:301,:,:))));
    Psioutput(yr,:)=Psires1;
    Psioutput2(yr,:)=Psires2;
end
%% Put it all together
% Make long vector of yearly max osfs
for i=1:length(Psioutput)-1
if i==1
    TS=[Psioutput(i,:),Psioutput(i+1,:)];
else
TS=[TS,Psioutput(i+1,:)];
end
end
for i=1:length(Psioutput2)-1
if i==1
TS2=[Psioutput2(i,:),Psioutput2(i+1,:)];
else
TS2=[TS2,Psioutput2(i+1,:)];
end
end
close all
%finish off the timeseries by adding in the 5 yeat tavs to the begining of
%the timeseries
%%PLOT TS
%Make a pretty plot and save it in script! 
x=[start:1:stop+4];%year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Residual Stream function','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Residual_osf_timeseries')
figure
plot(x,TS2,'k','linewidth',1.8);
title('Time series of Residual Stream function away from (Southern Boundary)','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Residual_osf_timeseries2')
end
if option ==2
    cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/Psi_dens/Rho_layers/
    for i=start:5:stop %5 year dumps  
    yr=(i-start)/5+1; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'Psi.nc'];
    [PsiTemp]=RSF3(fname,0); %Call my MOCcalc function no plot
    Psires1(1:5)=max(max(abs(PsiTemp(1:301,:,:))));
    Psires2(1:5)=max(max(abs(PsiTemp(50:301,:,:))));
    Psioutput(yr,:)=Psires1;
    Psioutput2(yr,:)=Psires2;
    end
%% Put it all together
% Make long vector of yearly max osfs
for i=1:length(Psioutput)-1
if i==1
    TS=[Psioutput(i,:),Psioutput(i+1,:)];
else
TS=[TS,Psioutput(i+1,:)];
end
end
for i=1:length(Psioutput2)-1
if i==1
TS2=[Psioutput2(i,:),Psioutput2(i+1,:)];
else
TS2=[TS2,Psioutput2(i+1,:)];
end
end
close all
%finish off the timeseries by adding in the 5 yeat tavs to the begining of
%the timeseries
%%PLOT TS
%Make a pretty plot and save it in script! 
x=[start:1:stop+4];%year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Residual Stream function','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Residual_osf_timeseries(dlayers)')
figure
plot(x,TS2,'k','linewidth',1.8);
title('Time series of Residual Stream function away from (Southern Boundary)','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-dpng','Residual_osf_timeseries2')
end
    
end