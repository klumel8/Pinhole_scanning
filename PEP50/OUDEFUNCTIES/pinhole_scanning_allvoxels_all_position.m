function [ in  ] = pinhole_scanning_allvoxels_all_position( pinholes, voxel_cor )
% This function takes all the informatien of the pinholes, such as location
% and orientation, togethr with thickness of the pinhole and the size of
% the angle of the pinhole. The second input are the coordinates of the
% different voxels.
% The output is a matrix booelan with the information if a pinhole sees a
% voxel.

%Hey
%creating cones section
%Matrix M (position vector cone, angle vector cone, maximum angle)
%for x y z phi theta radius alpha
M = pinholes;
C= zeros(size(M));
for i = 1:size(M,1)
C(i,1:3) = M(i,1:3) - (M(i,6)./tan(M(i,7))).*[cos(M(i,4)).*sin(M(i,5)) sin(M(i,4)).*sin(M(i,5)) cos(M(i,5))];
C(i,4:end) = M(i,4:end);
end 
%%  
N = size(voxel_cor,1);
%% 
%checking if point is inside cone(s)
int = ones(N,1);
in =ones(N,size(C,1));
for j = 1:size(C,1)
x0 = voxel_cor(:,1) - C(j,1);
y0 = voxel_cor(:,2) - C(j,2);
z0 = voxel_cor(:,3) - C(j,3);
v = [cos(C(j,4)).*sin(C(j,5)) sin(C(j,4)).*sin(C(j,5)) cos(C(j,5))]';
in(:,j) = (C(j,7) > acos( [x0 y0 z0]*v./(sqrt(sum(abs([x0 y0 z0]).^2,2))))); 
end 


end
