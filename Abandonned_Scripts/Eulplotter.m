function []=Eulplotter(fname,options)
%% A Plotter for Eulerian mean streamfunction 
%input file name and options
%%option 1 = area of interest plotted in standard plot
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/grid
lm=ncread('grid.nc','HFacS');
cd /noc/altix/scratch/hb1g13/MITgcm/nchannel/EulPsi/
if options == 1
Y=ncread(fname,'Yp1');
Eul=ncread(fname,'PsiVEL');
Eul(Eul==0)=NaN;
Eul=squeeze(nanmean(Eul));
Eul=squeeze(mean(Eul,3));
cd ../glued_state_files
Z=ncread('185-190all.nc','Z');
figure
contourf(Y(50:301)/1000,Z(1:15),Eul(50:301,1:15)'/10^6,15);
h=colorbar;
title('Eularian mean stream function','fontsize',12)
xlabel('Meridional distance S-N (km)','fontsize',12)
ylabel('Depth (m)','fontsize',12)
ylabel(h,'Transport (Sv)','fontsize',12)
print('MOCEul3','-dpng')
end
%%option 2 Overlay density for top 200m
if options ==2 
Y=ncread(fname,'Yp1');
Eul=ncread(fname,'PsiVEL');
%time average
Eultav=squeeze(mean(Eul,4));
%regrid to center x
for i=1:450
    Euln(i,:,:)=mean(Eultav(i:i+1,:,:));
end
Euln=Euln.*lm;
Euln(Euln==0)=NaN;
Euln=(squeeze(nanmean(Euln)));
cd ../glued_state_files
Z=ncread('185-190all.nc','Z');
figure
contourf(Y(1:301)/1000,Z(1:15),Euln(1:301,1:15)'/10^6,15);
h=colorbar;
title('Eularian mean stream function','fontsize',12)
xlabel('Meridional distance S-N (km)','fontsize',12)
ylabel('Depth (m)','fontsize',12)
ylabel(h,'Transport (Sv)','fontsize',12)
print('MOCEul3','-dpng')
end
if options > 2 || options < 1
    'Not an option....'
end

end