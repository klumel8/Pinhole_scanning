function quality = Main(gridrow, gridcol, translatedScanning, voxelCoords, part)



if true
    gridcir = CircleGrid4(gridrow,gridcol,part);
else
    load('circles.mat');    
end


quality = zeros(size(voxelCoords,1),1);

disp('Entering the voxelIntersectLoop')
quality = Intersections(gridrow,gridcol,voxelCoords,gridcir,translatedScanning);

end