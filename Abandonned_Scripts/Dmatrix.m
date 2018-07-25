%% To be added to Density_matrix.m 
clear
%Load Variables
%Change to Fname and put into function after testing
Y=ncread('220-225Psi.nc','Yp1'); %Y in km
Dens=csvread('dens'); % csv density bins
dz=ncread('220-225Psi.nc','LaHs1RHO'); %Layer thickness
Dwz=cumsum(dz,3); % Add the layer thickness to find the depth of density surface
Dens=Dens(2:49); % sill csv read starts with a 0 randomly, bin that and the diagnostics packages limits to 48 layers
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
save('Density_matrix2.mat','Dmatrix');
%% Plotting it
%D=squeeze(nanmean(Dmatrix)); %NB NANMEAN!!
%D=squeeze(mean(D,3));
%imagesc(Y,1:10:3500,D')