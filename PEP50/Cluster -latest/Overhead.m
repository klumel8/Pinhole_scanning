clear all

Lx = 1 
Ly = 30
Lz = 100
dx = 1
dy = 1
dz = 1
R_transaxial = 7
N_Pos_rev = 4.5
N_Pos_tot = 25

gridrow = 2001%input('Enter the gridrow'); %501
gridcol = 4000%input('Enter the gridcol'); %1000
division = 160%input('Enter the division');

parpool('local',16);


%{
gridrow
gridcol
division

parpool('local', 10)
%}

%{
gridrow = 1001;
gridcol = 2000;

Lx = 1; Ly = 1; Lz = 1;
dx = 1; dy = 1; dz = 1;
R_transaxial = 1;
N_Pos_rev = 20;
N_Pos_tot = 80;
%}

load('pinholes.mat');
scanning_input = [x y z phi -theta 0.5*d alpha];

%continuBeweging = false;   

translatedScanning = all_pinholepositions([Lx,Ly,Lz], R_transaxial,N_Pos_rev,N_Pos_tot,scanning_input);
%quiver3(translatedScanning(:,1),translatedScanning(:,2),translatedScanning(:,3),sin(translatedScanning(:,5)).*cos(translatedScanning(:,4)),sin(translatedScanning(:,5)).*sin(translatedScanning(:,4)),cos(translatedScanning(:,5)))
%translatedScanning = scanning_input;

voxelCoords = Voxel_coordinates(Lx,Ly,Lz,dx,dy,dz);
%voxelCoords = [0 0 0];

time = zeros(1,division);

partQuality = zeros(size(voxelCoords,1),division);
parfor go = 1:division
    %tic21
    %disp(go/division)
    disp('Entering Next division interation')
    part = [go division];
    partQuality(:,go) = Main(gridrow, gridcol, translatedScanning, voxelCoords, part);
    %partQuality(:,go);
    %time(1,go) = toc;
    %avg = mean(time(time~=0));
    %ETA_Overhead = (size(voxelCoords,1) - go)*avg/60;
    %fprintf("ETA_Overhead: %f",ETA_Overhead);
    
end
quality = mean(partQuality,2);
quality = reshape(quality,[Lz/dz,Ly/dy])';
save('savedQuality4','quality')
%colormap('jet');
%imagesc(quality);