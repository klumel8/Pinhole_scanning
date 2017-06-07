clear all
close all
clc
%% spiraatlje 
nvoxX=26;
nvoxY=26;
nvoxZ=81;
voxelsize=1;
radius_mouse=12.5;
x=[0 5 10 15 20 25];
mask = makemask(nvoxX,nvoxY,nvoxZ,voxelsize,radius_mouse);

for i= 1:length(x);                                     %calculating max, min, std, and mean at each number of positions for U-SPECT and VECTor with ST
steps=x(i);    
load (['Vector_150x75_' num2str(steps)])
load (['Uspect_150x75_' num2str(steps)])

Vmuis=mask.*Vquality;
Umuis=mask.*Uquality;
Vnonzero=nonzeros(Vmuis);
Unonzero=nonzeros(Umuis);

Umin=min(Unonzero);         Umax=max(Unonzero);
Vmin=min(Vnonzero);         Vmax=max(Vnonzero);
Vstd(i)=std(Vnonzero);      Vmean(i)=mean(Vnonzero);
Ustd(i)=std(Unonzero);      Umean(i)=mean(Unonzero);
end
%% plane                                                %calculating min, std and mean of U-SPECT with MPT
load Uspect_150x75_plane
steps=24;
plane=mask.*PUquality;
Pnonzero=nonzeros(plane);
Pmean=mean(Pnonzero);
Pstd=std(Pnonzero);
c=min(Pnonzero);

imagesc(squeeze(PUquality(13,:,:)));                    %plot of MPT in yz-plane of the middle slice
axis('equal');axis('tight');
colormap('jet');caxis([0.7 1]);
xlabel('z-axis (mm)'); ylabel('y-axis (mm)'); title(['Plane 24 steps ' num2str(x(i))]); colorbar('EastOutside');

%% quality curve of mean values of VECTor and U-SPECT with std as errorbar
figure(1);
errorbar(x,Vmean,Vstd); hold on; errorbar(x,Umean,Ustd); errorbar(steps,Pmean,Pstd,'o');hold on; legend('Vector','Uspect','Plane'); 
xlabel('# Bedpositions'); ylabel('Mean Quality'); title(' Mean Quality curve'); axis([-1 26 -0.1 1.2])
