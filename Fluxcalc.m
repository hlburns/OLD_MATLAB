%% Calculating Salt and Heat Fluxes
% The following function was written in Matlab(R) to calculate salt and 
% heat fluxes broken down into time mean and time varying components 
% This function is ran in parallel for loops for each 1/10^o latband and 
% each year.
function [TotalHeatFlux_PW, TotalSaltFlux_Sv,MeanHeatFlux_PW, MeanSaltFlux_Sv, EddyHeatFlux_PW, EddySaltFlux_Sv, OverturningSaltFlux_Sv, OverturningHeatFlux_PW, GyreSaltFlux_Sv, GyreHeatFlux_PW]=Fluxcalc(lat,yr)
%% This function calculates Salt and Heat Transports!
%From the following variables:
%S(X,Z,t) = Salinity at chosen latband
%V(X,Z,t) = NS Velocity at chosen latband
%U(X,Z,t) = EW Velocity at chosen latband
%Stave(X,Z) = Salinity timeavergaed over 5 years at chosen latband
%Vtave(X,Z) = NS Velocity timeaveraged over the 5 years at chosen latband
%Lat= Chosen latband
%Z = Depth of 33 layers of varying thickness
%Lon = 10th degree longitude
%Depth = Bathymetry data
%Area = Are of each cell
%Landmask = The original landmask from HFaCC file (% of cell open)
%% Add loctations of Time varying and Time averaged data
addpath (['/home/ocean1/helenb/Allyears/year',num2str(yr),'/sections/'])%Time varying
addpath /home/ocean1/helenb/Allyears/timeavs%5 year time averaged values
addpath /home/ocean1/helenb/Saltfluxes
load Sref %Sectional average saliniy in each lat band
load(['lat',num2str(lat),'tav.mat'])
Sref=Srefall(lat);
load('Zp.mat'); % Depth of the base of each box
Lat=(0.1*lat)-9.05; %convert numbers 1:800 to a Latitude
load (['e',num2str(Lat),'.mat']);
%%
%% Create Area Matrix

Area=zeros(1200,33);
for j=1:33
    for i=1:1200
L=(2*pi*cosd(lat)*6371*0.1)/360;%convert degrees of longitude into km
dz=Zp(1:33)-Zp(2:34); %Thickness of each box
A = L*10^3*dz(j);
Area(i,j)=A;
    end
end
As=ones(1200,33);
LA=L*As;
W=Landmask.*LA; 
Area=Landmask.*Area; %Add in the percentage of box open
%NaN for edge effects of contients to make matrix of 1s and NaNs;
Landmask(Landmask>0)=1; 
Landmask(Landmask==0)=NaN;
%% Appy Landmask
%Solving the hole in Greenland problem for HFaCC in year5;
Ve=zeros(1200,33,61);
Sa=zeros(1200,33,61);
T=zeros(1200,33,61);
if yr==5 && lat>61.15;
   V(V==0)=NaN;
   nm=isnan(V);
   snan=S.*nm;
   S=S-snan;
   S(S==0)=NaN;
   Tnan=Temp.*nm;
   Temp=Temp-Tnan;
   Temp(Temp==0)=NaN;
   Vtave(Vtave==0)=NaN;
   nm=isnan(Vtave);
   stavenan=Stave.*nm;
   Stave=Stave-stavenan;
   Stave(Stave==0)=NaN;
   Ttavenan=Ttave.*nm;
   Ttave=Ttave-Ttavenan;
   Ttave(Ttave==0)=NaN;   
end

for t=1:61
    V1=Landmask.*V(:,:,t);
    Ve(:,:,t)=V1;
    S1=Landmask.*S(:,:,t);
    Sa(:,:,t)=S1;
    T1=Landmask.*Temp(:,:,t);
    T(:,:,t)=T1;
end

%% Remove The Med and Hudson Bay and Normalise salt!
% The Med is very salty and the Hudson Bay is very fresh and are separate
% from the rest of the North Atlantic Basin.

if Lat > 50.05;
    Stave(1:400,1,:)=NaN;
    Ttave(1:400,1,:)=NaN;
    Vtave(1:400,1,:)=NaN;
    Ve(1:400,:,:)=NaN;
    Sa(1:400,:,:)=NaN;
    T(1:400,:,:)=NaN;
    
end
if Lat>21.05 && Lat<45.05
    Stave(1000:1200,1,:)=NaN;
    Ttave(1000:1200,1,:)=NaN;
    Vtave(1000:1200,1,:)=NaN;
    Ve(1000:1200,:,:)=NaN;
    Sa(1000:1200,:,:)=NaN;
    T(1000:1200,:,:)=NaN;
    if Lat>30.05 && Lat<35.05
        Stave(950:1200,:)=NaN;
        Ttave(950:1200,1,:)=NaN;
        Vtave(950:1200,1,:)=NaN;
        Ve(950:1200,:,:)=NaN;
        Sa(950:1200,:,:)=NaN;
        T(950:1200,:,:)=NaN;
    end
end
%Normalise Salinity:
% $$ S=(Sal-S_o)/S_o $$
% Where $S_o$ is the section average salinity.
Sal=(Sa-Sref)/Sref;

%% Create Values required for final total and time mean transports
% $$ \overline{VS}=\overline{\overline{V}\overline{S}}+\overline{V'S'} $$
VS=Ve.*Sal;
VT=Ve.*T;
Vetave=squeeze(Vtave).*Landmask;
Ttave=squeeze(Ttave).*Landmask;
Satave=(squeeze(Stave)-Sref)/Sref.*Landmask;
VStave=nanmean(VS,3);
VTtave=nanmean(VT,3);
%Total Flux: Timeav sum of VT product.*Area over the section.
% $$ \sum\sum \overline{VT} \Delta x \Delta z $$
TotalHeatFlux=nansum(VTtave(:,1:33).*(Area(:,1:33)));
TotalHeatFlux_PW=1030*3985*nansum(TotalHeatFlux)/(10^15);%Divide by 10^6 to get into Sv!

%Total Flux: Timeav sum of VS product.*Area over the section.
% $$ \sum\sum \overline{VS} \Delta x \Delta z $$
TotalSaltFlux=nansum(VStave(:,1:33).*(Area(:,1:33)));
TotalSaltFlux_Sv=nansum(TotalSaltFlux)/(10^6);%Divide by 10^6 to get into Sv!

%Mean_Flux: sum of Product of time averaged V and S * Area over the section
% $$ \sum\sum \overline{V}\overline{S} \Delta x \Delta z $$
MeanSaltFlux=nansum(Vetave(:,1:33).*Satave(:,1:33).*Area(:,1:33));
MeanSaltFlux_Sv=nansum(MeanSaltFlux)/10^6;
% $$ \sum\sum \overline{V}\overline{T} \Delta x \Delta z $$
MeanHeatFlux=nansum(Vetave(:,1:33).*Ttave(:,1:33).*Area(:,1:33));
MeanHeatFlux_PW=1030*3985*nansum(MeanHeatFlux)/10^15;

%Eddy Flux: Time av sum of V'S' * Area over section
% $$ \sum\sum \overline{V'S'} \Delta x \Delta z $$
% $$ \sum\sum \overline{V'T'} \Delta x \Delta z $$
VSprime=zeros(1200,33,61);
VTprime=zeros(1200,33,61);
%Calculate V' and S' at each timestep:
for t=1:61
Vprimet=Ve(:,:,t)-Vetave;
Sprimet=((Sal(:,:,t)-Satave));
VSprime(:,:,t)=(Vprimet.*Sprimet);
Tprimet=((T(:,:,t)-Ttave));
VTprime(:,:,t)=(Vprimet.*Tprimet);
end

%Time average V'S'
VSprimetav=nanmean(VSprime,3);
Eddy1=(VSprimetav);
Eddy=nansum((Eddy1(:,1:33)).*Area(:,1:33));
EddySaltFlux_Sv=(nansum(Eddy)/10^6);
VTprimetav=nanmean(VTprime,3);
Eddy2=(VTprimetav);
Eddy3=nansum((Eddy2(:,1:33)).*Area(:,1:33));
EddyHeatFlux_PW=1030*3985*(nansum(Eddy3)/10^15);


%Decomposing mean:
%Mean=Overturning+Gyre 
% $$\overline{V}\overline{S} = 

%% Overturing = $<V><S>$
Vzone=nanmean(Vetave,1);
Szone=nanmean(Satave,1);
Tzone=nanmean(Ttave,1);
Width=nansum(W.*(~isnan(Vetave)),1);
OverturningSaltFlux_Sv=(nansum(10^3*(Width(1:33).*Vzone(1:33).*Szone(1:33).*dz(1:33)'))/10^6);
OverturningHeatFlux_PW=1030*3985*(nansum(10^3*(Width(1:33).*Vzone(1:33).*Tzone(1:33).*dz(1:33)'))/10^15);
%% Gyre FLux:

%S* and V*

Vstar=zeros(1200,33);
Sstar=zeros(1200,33);
Tstar=zeros(1200,33);
VzoneTstar=zeros(1200,33);
VstarTzone=zeros(1200,33);
Vzonestar=zeros(1200,33);
Vstarzone=zeros(1200,33);
for i=1:1200;
    Vstar(i,:)=Vetave(i,:)-Vzone;
    Sstar(i,:)=Satave(i,:)-Szone;
    Tstar(i,:)=Ttave(i,:)-Tzone;
    Vstarzone(i,:)=Vstar(i,:).*Szone;
    Vzonestar(i,:)=Vzone.*Sstar(i,:);
    VstarTzone(i,:)=Vstar(i,:).*Tzone;
    VzoneTstar(i,:)=Vzone.*Tstar(i,:);
end
VStar=Vstar.*Sstar;
Gyre1_Sv=(nansum(nansum(VStar(:,1:33).*Area(:,1:33))))/10^6;
Gyre2_Sv=(nansum(nansum(Vstarzone(:,1:33).*Area(:,1:33))))/10^6;
Gyre3_Sv=(nansum(nansum(Vzonestar(:,1:33).*Area(:,1:33))))/10^6;
GyreSaltFlux_Sv=Gyre1_Sv+Gyre2_Sv+Gyre3_Sv;
VTstar=Vstar.*Tstar;
Gyre1_PW=1030*3985*(nansum(nansum(VTstar(:,1:33).*Area(:,1:33))))/10^15;
Gyre2_PW=1030*3985*(nansum(nansum(VstarTzone(:,1:33).*Area(:,1:33))))/10^15;
Gyre3_PW=1030*3985*(nansum(nansum(VzoneTstar(:,1:33).*Area(:,1:33))))/10^15;
GyreHeatFlux_PW=Gyre1_PW+Gyre2_PW+Gyre3_PW;
end