function []=Density_matrix(fname,options,layers)
%% Make Density matrix
%This function will calculate a density matrix from the layers output
%The inouts must be
%fname = String of the file name variable
%Option = 1 = Temp layers 2 = Density layers
%Layers = String containing the layers file
if options== 1
%Load Variables
%Change to Fname and put into function after testing
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/Psi_dens/
Y=ncread(fname,'Yp1'); %Y in km
Dens=csvread(layers); % csv density bins
dz=ncread(fname,'LaHs1TH'); %Layer thickness
dz(dz==0)=NaN;
Dwz=nancumsum(dz,3); % Add the layer thickness to find the depth of density surface
Dens=Dens(2:49);% silly csv read starts with a 0 randomly, bin that and the diagnostics packages limits to 48 layers
Dmatrix=zeros(450,301,140,4); 
%if matlabpool('size')==0
%    matlabpool open 4
%end
for t=1:4
    tic
    for i=1:450;
        for j=1:301
            for k=1:140
                %find the density in 10m intervals and set up a new matrix
                %of density with depth as the Z co-ordinate 
                A=Dens(find(Dwz(i,j,:,t)>25*k-25 & Dwz(i,j,:,t)<25*k));%contary to the hint here the find is neccessary
                if isempty(A)==1 %If there's no density values found at that depth then set to NaN
                   A=NaN; 
                elseif length(A)>1 %If there are many density vaules found that depth use the first occuring desity value
                    B=A(A~=0);
                    A=B(1);
                end
            Dmatrix(i,j,k,t)=A;
            end
        end 
    end
    toc
end
cd ~/PSIres_test_area/
save('Density_matrixT.mat','Dmatrix');
%% Plotting it
D=squeeze(nanmean(Dmatrix)); %NB NANMEAN!!
D=squeeze(mean(D,3));
imagesc(Y/1000,1:10:3500,D')
title('Density Field From Layers package','fontsize',12)
h=colorbar;
ylabel(h,'Density (kg/m^3)','fontsize',12)
xlabel('Meridional distance (km)','fontsize',12)
ylabel('Depth (m)','fontsize',12)
print('Density_field','-djpg')
end
if options == 2
    %Load Variables
%Change to Fname and put into function after testing
Y=ncread(fname,'Yp1'); %Y in km
Dens=csvread(layers); % csv density bins
dz=ncread(fname,'LaHs1RHO'); %Layer thickness
dz(dz==0)=NaN;
Dwz=nancumsum(dz,3); % Add the layer thickness to find the depth of density surface
Dens=Dens(2:49); % sill csv read starts with a 0 randomly, bin that and the diagnostics packages limits to 48 layers
RhoNil=1000;
tAlpha=2*10^-4;
tref=20.0;
Rho=RhoNil*(1-(tAlpha*(Temp-tref)));
Dmatrix=zeros(450,301,140,4); 
%if matlabpool('size')==0
%    matlabpool open 4
%end
for t=1:4
    tic
    for i=1:450;
        for j=1:301
            for k=1:140
                %find the density in 10m intervals and set up a new matrix
                %of density with depth as the Z co-ordinate 
                A=Rho(find(Dwz(i,j,:,t)>25*k-25 & Dwz(i,j,:,t)<25*k));%contary to the hint here the find is neccessary
                if isempty(A)==1 %If there's no density values found at that depth then set to NaN
                   A=NaN; 
                elseif length(A)>1 %If there are many density vaules found that depth use the first occuring desity value
                    B=A(A~=0);
                    A=B(1);
                end
            Dmatrix(i,j,k,t)=A;
            end
        end 
    end
    toc
end
cd ~/PSIres_test_area/
save('Density_matrixRho.mat','Dmatrix');
%% Plotting it
D=squeeze(nanmean(Dmatrix)); %NB NANMEAN!!
D=squeeze(mean(D,3));
imagesc(Y/1000,1:10:3500,D')
title('Density Field From Layers package','fontsize',12)
h=colorbar;
ylabel(h,'Density (kg/m^3)','fontsize',12)
xlabel('Meridional distance (km)','fontsize',12)
ylabel('Depth (m)','fontsize',12)
print('Density_field','-djpg')
end
if options > 2 && options <1
    'Option 1 = Density Layers'
    'Option 2 = Temp Layers'
end
   
end