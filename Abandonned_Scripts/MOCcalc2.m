function [ Psimax ] = MOCcalc2(Y,T,Psi)
% This will calculate the maximum overturning stream function. In depth coordinates
%
%Convert to Sv
Psi2=Psi/10^6;
%Find the max absolute value
Psimax=squeeze(max(max(max(abs(Psi2)))));
end  
