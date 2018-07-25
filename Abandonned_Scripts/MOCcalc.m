function [ Psimax ] = MOCcalc(Y,T,Psi)
% This will calculate the maximum overturning stream function. In depth coordinates
%
% This function requires latitude Y, Timesteps T, and the output stream
% function from the appropriate .nc file
%PsiVEL varies in X,Y,Z,T in m^3/s
%Convert to Sv
Psi2=Psi/10^6;
%Take zonal means - using the assumption the overturning varies much more
%N-S than it does E-W
Psizonal=squeeze(mean(Psi2));
%At every lat and time step find the max overturning stream function 
%I've bin the very surface and bottom values to remove ekman or strong
%bottom flow influences
if mod(length(Y),2) ~=0
   k=(length(Y)-1)/2;
else
   k=length(Y)/2;
end
Psimax1=zeros(k,length(T));
for i=1:length(T);
     for j=1:k; 
        Psimax1(j,i)=max(abs(Psizonal(k+j,1:15,i))); % Find the max transport.
     end
end
%No I have Psi varying with latitude and time
%I'm not sure where the maximum overturning is to be expected so I'll just
%find the max value at any latitude.
Psimax=(max(Psimax1));
clearvars -EXCEPT Psimax Y X Z T Psi
%Psimax now just varies with time - depending on timeaveraging used.
%Next step would be to concatenate these to give a time series.
end  
