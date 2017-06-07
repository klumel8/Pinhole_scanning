close all 
clear all
%% This script plots the middle slice of VECTor and U-SPECT for different numbers of bed positions in the yz-plane in one subplot figure
x=[0 5 10 15 20 25];            %vector of number of bed positions
min=0.6;                        %minimal value on colourbar
for i= 1:length(x);
steps=x(i);    
load (['Vector_150x75_' num2str(steps)])
load (['Uspect_150x75_' num2str(steps)])

subplot(length(x),2,(2*i-1));
imagesc(squeeze(Vquality(13,:,:)));
axis('equal');axis('tight');
colormap('jet');caxis([min 1]);
xlabel('z-axis (mm)'); ylabel('y-axis (mm)'); title(['VECTor ' num2str(x(i))]);
colorbar('EastOutside');

subplot(length(x),2,2*i);
imagesc(squeeze(Uquality(13,:,:)));
axis('equal');axis('tight');
colormap('jet');caxis([min 1]);
xlabel('z-axis (mm)'); ylabel('y-axis (mm)'); title(['U-SPECT ' num2str(x(i))]);
colorbar('EastOutside');
end