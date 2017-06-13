clear all

gridrow = 129;
gridcol = 256;

load('pinholes.mat');
theta = -theta;
scanning_input = [x y z phi theta 0.5*d alpha];
positions = [0 0 0];
voxelCoor = [0 0 0];

continuBeweging = false;

translatedScanning = translation(scanning_input, positions, gridrow, gridcol,continuBeweging);

if false
    gridcirTest = GridGenFloSas(gridrow,gridcol);
else
    load('circles.mat');
    
end
tic
%for itr = 1:150
singleVoxelPinholes = pinhole_scanning_1voxel(translatedScanning, voxelCoor);
quality = Intersections(gridrow,gridcol,singleVoxelPinholes,gridcirTest);
%end
toc