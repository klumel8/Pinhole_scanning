function quality = Main(gridrow, gridcol, scanning_input, scannerPositions, continuBeweging, voxelCoords)

gridrow = 401;
gridcol = 800;

load('pinholes.mat');
scanning_input = [x y z phi -theta 0.5*d alpha];
scannerPositions = [0 0 0; 0.1 0 0];

for i = 1:numel(voxelCoords)
    voxelCoor = voxelCoords(i)
end
voxelCoor = [0 0 0];

continuBeweging = false;

translatedScanning = translation(scanning_input, scannerPositions, gridrow, gridcol,continuBeweging);

if true
    gridcir = CircleGrid2(gridrow,gridcol);
else
    load('circles.mat');
    
end
tic
%for itr = 1:150
singleVoxelPinholes = pinhole_scanning_1voxel(translatedScanning, voxelCoor);
quality = Intersections(gridrow,gridcol,singleVoxelPinholes,gridcir);
%end
toc