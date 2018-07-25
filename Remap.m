%% Remap to Depth 
% Calls linearEOS.m and Resplotter.m (which calls RSF.m)
function [Psimap]=Remap(start,t,option)
%% This functions remaps Streamfunctions in (y,t) or (y,rho) into (y,Rho(z))
%Inputs:
%Fname = String of years e.g.225-230
%Option: 1= Temp layers 2=Density layers
%Calls RSF2 RSF3, linearEOSDens
close all
if nargin ~= 3 ;
   disp('input error');
   fprintf('input 3 variables: \n Start year \n time step \n option')
end
global Yc
global Zc
global T
global lmc
%% Density Layers
Comp='Iridis4';
if option ==1 
        OP='TEST';
        Comp='Nautilus';
elseif option==2
         OP='Sponge1';
elseif option==3
         OP='Sponge2';
elseif option==4
         OP='TEST';
end
fprintf(['option ',OP,' selected'])
cd(['/noc/msm/scratch/students/hb1g13/',Comp,'/',OP,'/Psi_dens/'])
%% Get T layers from linear EOS
Dens=csvread('Dens'); % CVS file used in data.layers
Rholayers=Dens(1:end-1); %Bin the last values
%% Get residual Streamfunction
stop=start; %Remap only one portion at a time
[~,PsiRho]=Resplottermsm(start,stop,5,option); %Call Residual Stream function calculator
close all %This function currently outputs plots - close them!
%% Get density field as F(y,z)
cd(['/noc/msm/scratch/students/hb1g13/',Comp,'/',OP,'/glued_state_files/'])
fname=[num2str(start),'-',num2str(stop+t),'all.nc'];
T=ncread(fname,'Temp');
cd(['/noc/msm/scratch/students/hb1g13/',Comp,'/',OP,'/grid'])
lmc=ncread('grid.nc','HFacC');
Zc=ncread('grid.nc','Z');
Yc=ncread('grid.nc','Y');
tref=csvread('Tref.out');
Rho=linearEOSDens(1,0,tref); % Calls linear EOS function asking for a density section and no plot
close all
%% Regridd Psi temp
%PsiRho(PsiRho==0)=NaN;
% Psi(x,yp1,zmd) to Psi(x,y,zmd)
PsiRho=squeeze(nanmean(PsiRho,3)); % time average
Psi=zeros(length(Yc),length(PsiRho(1,:)));
for i=1:length(Yc)
    Psi(i,:)=squeeze(nanmean(PsiRho(i:i+1,:),1));
end
%% Remap!
%I.E put Psi(x,y,zmd) to Psi(x,y,z);
%From Rho(y,z) 
%At each point find the density and find Psi for that latitude and density
%From Psi(y,Rho).
Psimap=zeros(length(Yc),length(Zc));
   for j=1:length(Yc)
        for k=1:length(Zc)
            D=Rho(j,k);
            P=Psi(j,:);
            SF=P(find(Rholayers>D-0.025 & Rholayers < D+0.025));
            if isempty(SF)==1 % The Density steps are bigger at depth
                  SF=P(find(Rholayers>=D-0.05 & Rholayers <=D+0.5));
            end

            if isempty(SF)==1 % The Density steps are bigger at depth
                  SF=P(find(Rholayers>=D-0.1 & Rholayers <= D+0.1));
            end
            if isempty(SF)==1 % The Density steps are bigger at depth
                  SF=P(find(Rholayers>=D-0.5 & Rholayers <= D+0.5));
            end
            if isempty(SF)==1
                    SF=NaN;
            end
            Psimap(j,k)=SF(1);
        end
   end
%% PLOT    
figure
contourf(Yc'/1000,Zc,Psimap',15);
cmax=max(max((Psimap)));
cmin=min(min((Psimap)));
colormap(b2r(cmin,cmax))
shading flat
title(['RMOC Remapped ',OP])
ylabel('Depth (m)')
xlabel('Meridional Distance (km)')
h=colorbar;
ylabel(h,'Transport (Sv)')
cd(['~/Figures/',Comp,'/',OP])
print([start,'yr_',OP,'-dpng'])
end
