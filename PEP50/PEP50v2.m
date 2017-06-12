% versie 1.0 
% PEP 50 -  pinhole quality check
% names:

clear all
load('Pinholes (1).mat')

% The area checked is a cube with a side 2L
L_area = 10;

% A voxel is a small cube with a side of length L
L_voxel = 1;

number_voxels = (L_area / L_voxel)^3;

% voxel coördinates 
[x_voxel,y_voxel,z_voxel] = Voxel_coordinates(L_area, L_voxel);

% fit the input information in an easy matrix
pinholes = [x,y,z,phi,theta,0.5*d, alpha] ;
more_positions
%returns a matrix for every voxel which pinhole it sees 
voxel_pinhole = pinhole_scanning(pinholes, [x_voxel',y_voxel',z_voxel']);

gridrow = 11;
gridcol = 10;
% Generates the circles on the sphere and puts the coordinates in a grid
circleGrid = CircleGridSasha(gridrow,gridcol);




