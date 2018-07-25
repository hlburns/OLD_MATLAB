function [respsi,maxofs]=res_psi(Fname)
%% Residual Stream Function
%read data
Y=ncread(Fname,'Yp1');
Z=ncread(Fname,'diag_levels');
T=ncread(Fname,'T');
V=ncread(Fname,'LaVa1RHO');
dz=ncread(Fname,'LaHs1RHO');
%load Dense2
%B=999.8:0.05:1000.2;
%C=1000.2:0.005:1000.31;
dx=6666.666; % 6.66km resolution
%for each time step
V(V==0)=NaN;
dz(dz==0)=NaN;
psi=zeros(301,length(Z));
psi_r=zeros(301,length(Z),5);
for t=1:length(T);
    for y=1:length(Y); % for each lat
        v_lat=squeeze(V(:,y,:,t));
        dz_lat=squeeze(dz(:,y,:,t));
        Vdz=v_lat.*dz_lat;
        %intVdz=cumsum(Vdz,2);
        intintVdzdx=nansum(dx*Vdz);
        psi(y,:)=intintVdzdx;
    end
    psi_r(:,:,t)=psi;
end
respsi=mean(psi_r,3);
%Plotting
%uimagesc(Y,B,respsi(:,2:20)/10^6)
maxofs=max(max(respsi))/10^6;
end