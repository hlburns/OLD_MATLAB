Psi=zeros(301,24,450); 
for j=1:301
     V(V==0)=NaN;
     V_lat=-V(:,1,:);
     v_lat = squeeze(V_lat); % V(Lon,depth,time)
     v_lat = permute(v_lat,[2 1]); % Rearrage V to V(Depth, Long)
     dx=6666.666; % dx is the width of the boxes at each individual latitude.
     dz=Z(1:23)-Z(2:24);
     dz=[0-Z(1);dz];
     osf=ones(length(Z),length(X)); % Pre-allocating a matrix of overturning stream function.
     osf(end,:)=0; % Forces osf to 0.
     for i=length(dz)-1:-1:1 % Start at the bottom and work up. 
         osf(i,:)=(nansum(squeeze(v_lat(i,:)).*(dx*ones(1,size(v_lat,2))).*dz(i),2))'+osf(i+1,:); % Work out the transports for each depth. 
     end
     Psi(j,:,:)=osf(:,:); % Bin the bottom variable 
        %osf_max=max(osf(find(-Zp>500 & -Zp<2500),:),[],1); % Find the max transport.
        %MOCval(t) = mean(osf_max); %This is the way to get it according to Bugnion et al.
        %clearvars -EXCEPT MOCval t Zp Vel Landmask v_lat V1 V osf
end