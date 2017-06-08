% versie 1.0 
% PEP 50 -  pinhole quality check
% names:


% The area checked is a cube with a side 2L
L_area = 1;

% A voxel is a small cube with a side of length L
L_voxel = 0.1;

number_voxels = (L_area / L_voxel)^3;

% voxel coördinates 
[x_voxel,y_voxel,z_voxel] = Voxel_coordinates(L_area, L_voxel);
