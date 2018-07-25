function [Psi]=RSF(Fname)
%% Residual Stream function calculator
V=ncread(Fname,'LaVa1RHO');% Velocity at density points
dz=ncread(Fname,'LaHs1RHO');% Thickes of deisty layer
%First integrate Vdz from depth to the surface
V(V==0)=NaN;
dz(dz==0)=NaN;
Vdz=V.*dz;
for i=1:48
intVdz(:,:,i,:)=nansum(Vdz(:,:,1:i,:),3);
end
%Then integrate Int (IntVdz)dx
%dx = 6.66km res
Psi=squeeze(nanmean(squeeze(nanmean(6666.666*(intVdz))),3));


end