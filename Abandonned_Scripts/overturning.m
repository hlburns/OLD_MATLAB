function [Psi] = overturning(GRID,v,varargin);
%psi=overturning(GRID,v);
%psi=overturning(GRID,v,mask);
%
%Calculates the volume overturning stream function (m^3/s).
%e.g.
% >> GRID=loadgrid;
% >> STATE=loadstate;
% >> psi=overturning(GRID,STATE.V,GRID.matl);
% >> contourf(GRID.yg(1,:),GRID.rf,psi'/1e6);colorbar
%
%Written by adcroft@mit.edu, 2002
%$Header:

if nargin==3
 mapmsk=varargin{:};
elseif nargin==2
 mapmsk=max( GRID.mskc(:,:,1), GRID.mskc(:,:,end) );
else
 error('Wrong number of arguments')
end
mapmsk=mapmsk.*mapmsk(:,[1 1:end-1]);

N=size(GRID.hfacs);
nxy=prod(size(GRID.rac));
nr=prod(size(GRID.drf));

DRF=spdiags(GRID.drf',0,nr,nr);
dz=reshape(GRID.hfacs,[nxy nr])*DRF;
area=dz.*( (GRID.dxg(:).*mapmsk(:))*ones(1,nr) );
area=reshape(area, N);

V=squeeze( sum(v.*area,1) );

Psi=zeros(N(2),nr+1);
for k=nr:-1:1;
 Psi(:,k)=Psi(:,k+1)-V(:,k);
end

vmsk=squeeze( sum(GRID.hfacs.*area,1) );
vmsk( find(vmsk~=0) )=1;
pmsk=vmsk(:,[1 1:end]).*vmsk(:,[1:end end]);

pmsk( find(pmsk==0) )=NaN;
Psi=Psi.*pmsk(:,[1 1:end-1]);
