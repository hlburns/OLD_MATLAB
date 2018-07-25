%% MOC calculator
%clear
%close all
load('26.5Nvelmonthlymeans.mat') % Load 26.6N
load 'Zp.mat' % Depth of the base of each box.
Vel=permute(Vel,[2 3 1]);% Make velocitites positive south and adjust for % of grid box open
Landmask=permute(Landmask,[2 1]);
for t = 1:305; % For each time step
        %% Calculate MOC strength at 26.666N. 
        V_lat=-Vel(:,:,t).*Landmask;
        lat_1 = 26.55; % Take latitudes of interest. 
        v_lat = V_lat(:,:); % V(Lon,depth,time)
        %v_lat = permute(v_lat,[2 1]); % Rearrage V to V(Depth, Long)
        dx=111320.*cosd(lat_1)./10; % dx is the width of the boxes at each individual latitude.
        dz=Zp(1:33)-Zp(2:34); % dz is a list of box thicknesses.  

        osf=ones(length(Zp)+1,length(lat_1)); % Pre-allocating a matrix of overturning stream function.
        osf(end,:)=0; % Forces osf to 0.

        for i=length(dz):-1:1 % Start at the bottom and work up. 
            osf(i,:)=(sum(squeeze(v_lat(i,:)).*(dx*ones(1,size(v_lat,2))).*dz(i),2))'+osf(i+1,:); % Work out the transports for each depth. 
        end
        osf=osf(1:end-1,:); % Bin the bottom variable 
        osf_max=max(osf(find(-Zp>500 & -Zp<2500),:),[],1); % Find the max transport.
        MOCval(t) = mean(osf_max); %This is the way to get it according to Bugnion et al.
        clearvars -EXCEPT MOCval t Zp Vel Landmask v_lat V1 V osf
    
end