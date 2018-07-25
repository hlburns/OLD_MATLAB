function [TS]=EulTS(Nyrs1,Nyrs5,start,stop)
%% Eularian Stream function timeseries
%150-170years =5 year average
%170-200years = yearly average
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/EulPsi/
%Define outputmatrix size (no. years,5 yearly avs per file) 
Psioutput=zeros(Nyrs1/5,5); 
Psioutput2=zeros(Nyrs5/5,1);
%% Read .nc files and put into MOCcalc
for i=start:5:stop %5 year dumps  
    yr=(i/10)*2-29; %incremental steps from 1 for filling in matix
    yr2=i+5; %5 year steps to make file name
    %generate filename 
    fname=[num2str(i),'-',num2str(yr2),'EulPsi.nc'];
    T=ncread(fname,'T');
    %Some times there are 5 year averages rather than yearly averages
    if length(T)==1 
    Y=ncread(fname,'Yp1');
    Psi=ncread(fname,'PsiVEL');
    PsiEul=MOCcalc(Y,T,Psi); %Call my MOCcalc function
    Psioutput2(yr,:)=PsiEul;
    continue
    end
    yr3=yr-length(Psioutput2); %restart to 1 when yearly tavs start
    Y=ncread(fname,'Yp1');
    Psi=ncread(fname,'PsiVEL');
    PsiEul=MOCcalc(Y,T,Psi); %Call my MOCcalc function
    Psioutput(yr3,:)=PsiEul;
end
%% Put it all together
% Make long vector of yearly max osfs
for i=1:length(Psioutput)-1
if i==1
A=[Psioutput(i,:),Psioutput(i+1,:)];
else
A=[A,Psioutput(i+1,:)];
end
end
%finish off the timeseries by adding in the 5 yeat tavs to the begining of
%the timeseries
TS=[Psioutput2',A];
%%PLOT TS
%Make a pretty plot and save it in script! 
x=[start+1:5:start+Nyrs5-1,start+Nyrs5:1:stop+5-1]; %year axis in 5 year time steps and then yearly steps
plot(x,TS,'k','linewidth',1.8);
title('Time series of Eularian Stream function','fontsize',12)
xlabel('Years','fontsize',12)
ylabel('Streamfunction (Sv)','fontsize',12)
cd /noc/users/hb1g13/Figures/Nautilus/Run1
print('-deps','Eulerian_Mean_timeseries2')
end
 
 