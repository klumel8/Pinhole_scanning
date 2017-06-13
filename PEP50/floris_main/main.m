%versie 1.0
grid_row = 129;
grid_col = 256;

load('pinholes.mat')
theta = -theta;
c_input = [x y z phi theta 0.5*d alpha];
v_xyz = rand(1,3)*10 - 5;
%v_xyz = [0 0 0];
hit_pin = pinhole_scanning_1voxel(c_input, v_xyz);
grid_covered = circle_area(hit_pin, grid_row, grid_col);
%grid_circles = CircleGrid(grid_row, grid_col);
load('GC.mat');
tic
Intersect(GC, grid_covered);
toc
