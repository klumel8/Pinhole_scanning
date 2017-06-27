clear all

Lx = 1; Ly = 25; Lz = 80 %define the size(mm) of the area youre looking at (around the origin).

dx = 1; dy = 1; dz = 1;  %define the size(mm) of each voxel in the x,y,z, directions.

R_transaxial = 7
N_Pos_rev = 4.5
N_Pos_tot = 15

%give the grid resolution as input.
gridrow = input('Enter the gridrow');
gridcol = input('Enter the gridcol');

%Make divisions to split up the calculations and allow for bigger calculations.
division = input('Enter the division');

%load all the pinhole data, returns: x, y, z, phi, theta, d, alpha.
load('pinholes.mat');

%A dummy variable that will shorten the code later on.
scanning_input = [x y z phi -theta 0.5*d alpha];

%Add the translated coordinates to the pinhole coordinates.
translatedScanning = all_pinholepositions([Lx,Ly,Lz], R_transaxial,N_Pos_rev,N_Pos_tot,scanning_input);

%VoxelCoords is a matrix of containing the [x y z] of each voxel inside the area.
voxelCoords = Voxel_coordinates(Lx,Ly,Lz,dx,dy,dz);

%Pre-Allocate an array for the ETA function.
time = zeros(1,division);

%Pre-Allocate a matrix to determine the quality later on.
partQuality = zeros(size(voxelCoords,1),division);

%Loop through all calculations.
for go = 1:division
    %start the tic for the ETA function
    if ETA
        tic
    end
    
    %Display go/division, this gives a rough estimation of the progress.
    disp(go/division);
    disp('Entering Next division interation');
    
    %part is the current division of the total division.
    part = [go division];
    
    %calculate the quality for the division.
    partQuality(:,go) = Main(gridrow, gridcol, translatedScanning, voxelCoords, part);
    
    %The estimated time until completion function (minutes
    if ETA
        time(1,go) = toc; avg = mean(time(time~=0)); ETA_Overhead = (size(voxelCoords,1) - go)*avg/60;
        fprintf("ETA_Overhead: %f minutes and %f seconds",floor(ETA_Overhead),round((ETA_Overhead-floor(ETA_Overhead))*60));
    end
    
end
quality = mean(partQuality,2);
quality = reshape(quality,[Lz/dz,Ly/dy])'
