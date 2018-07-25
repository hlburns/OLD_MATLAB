%% SPONGE

%% DOMAIN, Depth, Q, tau, flat bottom = Abernathey 2011

% This is a matlab script that generates the input data
% This version is designed for a Southern Ocean simulation
% Previously commented binary outputs reinstated,
% in case program has problems with netCDF input.
% All binary output now in 'native' format, rather
% than explicitly big-endian = 'b' as originally.

clear all
close all
addpath '/noc/users/jeff/packages/mexcdf/mexcdf.r4040/mexnc'
addpath '/noc/users/jeff/packages/mexcdf/mexcdf.r4040/snctools'
% Dimensions of grid
nx=200;
ny=400;
nz=30;

% Nominal depth of model (meters)
H=2985;

% Size of domain in y & x
Ly=2000e3;
Lx=1000e3;

% Horizontal resolution (m) in y first
dymax = 5000.00 ; 
dy = ones(ny,1)*dymax;

%% Make x resolution equal to max of y resolution
dx = ones(nx,1)*dymax;

%% Write these to file
fid=fopen('delY','w','native'); fwrite(fid,dy,'real*8'); fclose(fid);
fid=fopen('delX','w','native'); fwrite(fid,dx,'real*8'); fclose(fid);

%% Calculate the vector x with gridpoints in the dead centre.
x=cumsum(dx)-dx/2;
y=cumsum(dy)-dy/2;
x=x-Lx/2;
y=y-Ly/2;
y(end) + dymax/2;

[Y,X]=meshgrid(y,x);
 
% Surface heat flux (W/m^2) -- three different realisations
Qo=10;
% Q2 = zeros(size(Y));
% y2=cumsum(dy);
% for i=1:nx
% Q2(i,1:334)=Qo*(cos(3*pi*y2(1:334)/(Ly)));
% end
Q2=Qo*sin(3*pi*Y/Ly);
Q2(:,334:end)=0;
fid=fopen('Qsurface.2','w','native'); fwrite(fid,Q2,'real*8'); fclose(fid);
% THIS BIT WRITES OUT THE VARIABLE IN NETCDF FORMAT
nc_create_empty('Qsurface.2.nc')
nc_add_dimension('Qsurface.2.nc','X',nx)
nc_add_dimension('Qsurface.2.nc','Y',ny)
varstruct.Name = 'Q';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Y','X'};
varstruct.Attribute = {};
nc_addvar('Qsurface.2.nc',varstruct)
nc_varput('Qsurface.2.nc','Q',Q2') 

%% Make topography -- flat to start with
h = -H*ones(nx,ny);
%% First, add a ridge-like feature.
%h = h+(2500 + 300*sin(10*pi*Y/Ly) + 400*sin(8*pi*Y/Ly)+ 300*sin(25*pi*Y/Ly) ).*sech((X-0.2*Y+3e5)/1.2e5);
%% Now add another ridge!!
%h = h+(2000 + 600*sin(11*pi*Y/Ly) + 300*sin(7*pi*Y/Ly)+ 500*sin(21*pi*Y/Ly) ).*sech((X+0.1*Y+1.5e6)/1.2e5);
%% Next add a continental slope next to Antarctica
%h = max(h,-0.01*Ly/2-0.01*Y);
%% Make sure channel is closed.
h(:,end) = 0;
h(:,1) = 0;
fid=fopen('topog','w','native'); fwrite(fid,h,'real*8'); fclose(fid);
% THIS BIT WRITES OUT THE VARIABLE IN NETCDF FORMAT
nc_create_empty('topog.nc')
nc_add_dimension('topog.nc','X',nx)
nc_add_dimension('topog.nc','Y',ny)
varstruct.Name = 'h';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Y','X'};
varstruct.Attribute = {};
nc_addvar('topog.nc',varstruct)
nc_varput('topog.nc','h',h') 

%% now do wind stress -- again, three possibilities
taumax = 0.2;
tau=taumax*max(cos(pi*Y/Ly),0);
%tau2=zeros(nx+1,ny);
%for i=1:nx+1
%    tau2(i,:)=tau;
%end
fid=fopen('windx.2','w','native'); fwrite(fid,tau,'real*8'); fclose(fid);
% THIS BIT WRITES OUT THE VARIABLE IN NETCDF FORMAT
tau2x=[tau; tau(1,:)];
nc_create_empty('wind.2.nc')
nc_add_dimension('wind.2.nc','X',nx)
nc_add_dimension('wind.2.nc','Xp1',nx+1)
nc_add_dimension('wind.2.nc','Y',ny)
nc_add_dimension('wind.2.nc','Yp1',ny+1)
varstruct.Name = 'taux';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Y','Xp1'};
varstruct.Attribute = {};
nc_addvar('wind.2.nc',varstruct)
nc_varput('wind.2.nc','taux',tau2x') 
tau2y = zeros(nx,ny+1);
varstruct.Name = 'tauy';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Yp1','X'};
varstruct.Attribute = {};
nc_addvar('wind.2.nc',varstruct)
nc_varput('wind.2.nc','tauy',tau2y') 

%% Now a new 3D diffusion file.
diffusi = zeros(nx,ny,nz);
for kk = 1:nz
  diffusi(:,:,kk) =  max(20*0.005/Ly*Y-0.045,0);
  diffnew(kk,:,:) =  squeeze(diffusi(:,:,kk))';  
end
fid=fopen('diffusi','w','native'); fwrite(fid,diffusi,'real*8'); fclose(fid);
% THIS BIT WRITES OUT THE VARIABLE IN NETCDF FORMAT
nc_create_empty('diffusi.nc')
nc_add_dimension('diffusi.nc','X',nx)
nc_add_dimension('diffusi.nc','Y',ny)
nc_add_dimension('diffusi.nc','Z',nz)
varstruct.Name = 'krho';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Z','Y','X'};
varstruct.Attribute = {};
nc_addvar('diffusi.nc',varstruct)
nc_varput('diffusi.nc','krho',diffnew) 

%afig(3)
figure(3)
subplot(221),plot(y/1e3,Q2(1,:))
xlabel('y (km)')
ylabel('Q (W/m^2)')
title('Surface heat flux (W/m^2)')
legend('Q1')
axis tight

subplot(222),plot(y/1e3,tau(1,:))
xlabel('y (km)')
ylabel('\tau (N/m^2)')
title('Wind Stress (N/m^2)')
legend('\tau 2')
axis tight

subplot(224),plot(y/1e3,h(1,:),y/1e3,h(ceil(nx/4),:),y/1e3,h(ceil(nx/3),:),y/1e3,h(ceil(nx/2),:))
xlabel('y (km)')
ylabel('h (m)')
axis tight
legend(sprintf('x = %5.0f km',x(1)/1e3),sprintf('x = %5.0f km',x(ceil(nx/4))/1e3),sprintf('x = %5.0f km',x(ceil(nx/3))/1e3),sprintf('x = %5.0f km',x(ceil(nx/2))/1e3))
title('Topography (m)')

%subplot(224),pcolor(x/1e3,y/1e3,h')
%shading flat
%colorbar
%set(gca,'dataaspectratio',[1 1 1])
%xlabel('x (km)') 
%ylabel('y (km)') 
%hold on
%contour(x/1e3,y/1e3,h',[0 0],'k')
%title('Topography (m)')

subplot(223),plot(y/1e3,squeeze(diffusi(1,:,1)),y/1e3,squeeze(diffusi(1,:,4)),y/1e3,squeeze(diffusi(1,:,7)),y/1e3,squeeze(diffusi(1,:,11)),y/1e3,squeeze(diffusi(1,:,12)))
xlabel('y (km)') 
ylabel('Diffusion (m^2/s)') 
title('Additional Vertical Diffusion (m^2/s)') 
set(gca,'ylim',[-0.001 0.006],'xlim',[-1250 1250])
legend('surface','300m','600m','1400m','deep',2)
%axis tight 

print -dpdf plot_Sponge_helen.pdf

%% Create Sponge stuff
%- rbcs mask & restoring tracer field:
msk=zeros(nx,ny,nz);
for i=1:nz
  for j=ny-20:ny
    msk(:,j,i)=(Y(:,j)-Y(:,ny-20))./(Y(:,ny)-Y(:,ny-20));
  end
end
%ieee='b';
prec='real*8';
fid=fopen('rbcs_mask.bin','w','ieee-le'); fwrite(fid,msk,prec); fclose(fid);
msk2=zeros(nx,ny,nz);
msk2(:,ny-20:ny,:)=1;
%ieee='b';
prec='real*8';
fid=fopen('rbcs_mask2.bin','w','ieee-le'); fwrite(fid,msk2,prec); fclose(fid);                                                                  
nc_create_empty('RBCS.nc')                                                                          
nc_add_dimension('RBCS.nc','X',nx)                                                                  
nc_add_dimension('RBCS.nc','Y',ny)                                                                  
nc_add_dimension('RBCS.nc','Z',nz)                                                                  
varstruct.Name = 'mask';                                                                             
varstruct.Nctype = 'double';                                                                         
varstruct.Dimension = {'Z','Y','X'};                                                                 
varstruct.Attribute = {};                                                                            
msk=permute(msk,[3 2 1]);
nc_addvar('RBCS.nc',varstruct)                                                                      
nc_varput('RBCS.nc','mask',msk) 
figure
plot(y/1e3,squeeze(msk(1,:,1)))
%Create temp_sponge.bin
% Define temperature profile from Abernathey 2011
% --------------------------
% Stratification (uniform)
N = 1.0e3; 
deltaT=8;
Tref = zeros(nz,1);
a=[5,22.5,60];
b=135:100:2435;
c=2535:150:2885;
z=-[a,b,c];
for k=1:nz
  Tref(k) = deltaT*(exp(z(k)/N)-exp(-H/N))/(1-exp(-H/N));
end
figure
plot(Tref,z)
t = zeros(nx,ny,nz);
t2=zeros(nz,ny,nx);
for k=1:nz
  t(:,:,k) = Tref(k);
  t2(k,:,:)= Tref(k);
end
fid = fopen('Tref.file','w');
fprintf(fid,'Tref = \n');
fprintf(fid,'  %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',Tref);
fclose(fid);
% Also write binary copy
fid = fopen('Tref.file.bin','w','ieee-le');
fwrite(fid,Tref,prec); fclose(fid);

nc_create_empty('Tinit.nc')
nc_add_dimension('Tinit.nc','X',nx)
nc_add_dimension('Tinit.nc','Y',ny)
nc_add_dimension('Tinit.nc','Z',nz)
varstruct.Name = 'temp';
varstruct.Nctype = 'double';
varstruct.Dimension = {'Z','Y','X'};
varstruct.Attribute = {};
nc_addvar('Tinit.nc',varstruct)
nc_varput('Tinit.nc','temp',t2) 
% Also write a formatted copy of Tref to a file,
% for use in producing the data file for runs
fid = fopen('Tref.out','w');
fprintf(fid,'Tref = \n');
fprintf(fid,'  %10.6f,%10.6f,%10.6f,%10.6f,%10.6f,\n',Tref);
fclose(fid);
% Also write binary copy
fid = fopen('Temp_sponge.bin','w','ieee-le');
fwrite(fid,t,prec); fclose(fid);






