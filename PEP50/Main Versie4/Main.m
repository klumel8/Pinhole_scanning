function quality = Main(gridrow, gridcol, translatedScanning, voxelCoords, part)



if true
    gridcir = CircleGrid4(gridrow,gridcol,part);
else
    load('circles.mat');    
end


quality = zeros(size(voxelCoords,1),1);

disp('Entering the voxelIntersectLoop')
for i = 1:size(voxelCoords,1)
    voxelCoor = voxelCoords(i,1:3);
singleVoxelPinholes = pinhole_scanning_1voxel(translatedScanning, voxelCoor);
quality(i) = Intersections(gridrow,gridcol,singleVoxelPinholes,gridcir);
end

end