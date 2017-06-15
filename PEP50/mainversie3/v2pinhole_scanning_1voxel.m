function [ hit  ] = pinhole_scanning_1voxel( pinholes, voxel_cor )
% This function takes all the informatien of the pinholes, such as location
% and orientation, together with thickness of the pinhole and the size of
% the angle of the pinhole. The second input are the coordinates of the
% voxel.
% The output is a matrix  with the information if a pinhole sees a
% voxel.

%for M =[x y z phi theta radius alpha]
%x,y,z are coordinates pinholes
%phi, theta are angles of normal vector of pinhole
%radius is the radius of the pinhole
%alpha is half of the angle of the pinhole
M = pinholes;
C= zeros(size(M));
for i = 1:size(M,1)
C(i,1:3) = M(i,1:3) - (M(i,6)./tan(M(i,7))).*[cos(M(i,4)).*sin(M(i,5)) sin(M(i,4)).*sin(M(i,5)) cos(M(i,5))];
C(i,4:end) = M(i,4:end);
end 
%% 
%checking if point is inside cone(s)
hit(:,1) = [voxel_cor,[0 0 0 0]];
p = 2;
for j = 1:size(C,1);
x0 = voxel_cor(1) - C(j,1);
y0 = voxel_cor(2) - C(j,2);
z0 = voxel_cor(3) - C(j,3);
v = [cos(C(j,4)).*sin(C(j,5)) sin(C(j,4)).*sin(C(j,5)) cos(C(j,5))]';
if (C(j,7) >= abs(acos( [x0 y0 z0]*v./(sqrt(sum(abs([x0 y0 z0]).^2,2)))))) ... 
    && sqrt(sum(abs([x0 y0 z0]).^2,2)) >= sqrt(sum(abs([M(j,1) - C(j,1) M(j,2) - C(j,2) M(j,3) - C(j,3)]).^2,2)) 
            hit(:,p) = M(j,:)';
            p = p + 1;       
end 
end 
end