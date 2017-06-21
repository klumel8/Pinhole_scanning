clear all

Lx = 10
Ly = 10
Lz = 10
dx = 1
dy = 1
dz = 1
R_transaxial = 3
N_Pos_rev = 8
N_Pos_tot = 8

gridrow = input('Enter the gridrow') %501
gridcol = input('Enter the gridcol') %1000
division = input('Enter the division')


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
%translatedScanning = scanning_input;

voxelCoords = Voxel_coordinates(Lx,Ly,Lz,dx,dy,dz);
%voxelCoords = [0 0 0];

partQuality = zeros(size(voxelCoords,1),division);
for go = 1:division
    disp(go/division)
    disp('Entering Next division interation')
    part = [go division];
    partQuality(:,go) = Main(gridrow, gridcol, translatedScanning, voxelCoords, part);
    partQuality(:,go)
end
quality = mean(partQuality,2);