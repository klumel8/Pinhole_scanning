function [ x,y,z ] = Voxel_coordinates( L_area, L_voxel )
% Returns the x,y,z coorinates of all the voxels as vectors
% The input arguments are the size of the area and the voxels

xt = linspace(-L_area,L_area, L_area/L_voxel);
yt = linspace(-L_area,L_area, L_area/L_voxel);
zt = linspace(-L_area,L_area, L_area/L_voxel);
t = combvec(xt,yt, zt);
x = t(1,:);
y = t(1,:);
z = t(1,:);
end

