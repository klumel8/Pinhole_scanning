%% Quiver box plot
close all;
Lz=80; Lx=25; Ly=25; dx=1; dy=1; dz=1;

figure(1);

hold on
quiver3(rpin(:,1),rpin(:,2),rpin(:,3),normpin(:,1),normpin(:,2),normpin(:,3));

size = [Lx Ly 20];
origin = [-Lx/2 -Ly/2 -20/2];
x=([0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]-0.5)*size(1)+origin(1)+size(1)/2;
y=([0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]-0.5)*size(2)+origin(2)+size(2)/2;
z=([0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]-0.5)*size(3)+origin(3)+size(3)/2;
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),'w');
    set(h,'edgecolor','k')
end 

axis image
xlabel('x-axis (mm)'); ylabel('y-axis (mm)'); zlabel('z-axis (mm)')

hold off

%% Image generation
pagex = 13;                 %the slice displayed in yz-plane
pagez = 40;                 %the slice displayed in xy-plane
a1 = -Lx/2:dx:Lx/2;         %size of axis
a2 = -Ly/2:dy:Ly/2;
min=0.8;                    %minimal value on colourbar

subplot(2,2,1);
imagesc(squeeze(Vquality(:,:,pagez)));
axis('tight');colormap('jet');caxis([min 1]);colorbar('EastOutside');
xlabel('x-axis (mm)'); ylabel('y-axis (mm)'); title('VECTor xy ');

subplot(2,2,3);
imagesc(squeeze(Vquality(pagex,:,:)));
axis('equal');axis('tight');colormap('jet');caxis([min 1]);
xlabel('z-axis (mm)'); ylabel('y-axis (mm)'); title('Vector yz');

subplot(2,2,2);
imagesc(squeeze(Uquality(:,:,pagez)));
xlabel('x-axis (mm)'); ylabel('y-axis (mm)'); title('Uspect xy');
axis('tight');colormap('jet');caxis([min 1]);colorbar('EastOutside');

subplot(2,2,4);
imagesc(squeeze(Uquality(pagex,:,:)));
axis('equal');axis('tight');colormap('jet');caxis([min 1]);
xlabel('z-axis (mm)'); ylabel('y-axis (mm)'); title('Uspect yz');

