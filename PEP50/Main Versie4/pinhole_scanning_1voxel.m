function [ in  ] = pinhole_scanning_1voxel( pinholes, voxel_cor )

% This function takes all the informatien of the pinholes, such as location
% and orientation, togethr with thickness of the pinhole and the size of
% the angle of the pinhole. The second input are the coordinates of the
% different voxels.
% The output is a matrix  with the information if a pinhole sees a
% voxel.

%Hey
%creating cones section
%Matrix M (position vector cone, angle vector cone, maximum angle)
%for x y z phi theta radius alpha

    M = pinholes;
    C= zeros(size(M));
    C(:,1:3) = M(:,1:3) - (M(:,6)./tan(M(:,7))).*[cos(M(:,4)).*sin(M(:,5)) sin(M(:,4)).*sin(M(:,5)) cos(M(:,5))];
    C(:,4:end) =M(:,4:end);
    %{
for i = 1:size(M,1)
C(i,1:3) = M(i,1:3) - (M(i,6)./tan(M(i,7))).*[cos(M(i,4)).*sin(M(i,5)) sin(M(i,4)).*sin(M(i,5)) cos(M(i,5))];
C(i,4:end) = M(i,4:end);
end
    %}
    %%
    
    %checking if point is inside cone(s)
    
    in(:,1) = [voxel_cor,[0 0 0 0]];
    p = 2;
    %{
    xtest = voxel_cor(1) - C(:,1);
    ytest = voxel_cor(2) - C(:,2);
    ztest = voxel_cor(3) - C(:,3);
    vtest = [cos(C(:,4)).*sin(C(:,5)) sin(C(:,4)).*sin(C(:,5)) cos(C(:,5))]';
    chi = C(:,7)> acos( [xtest ytest ztest]*vtest./(sqrt(sum(abs([xtest ytest ztest]).^2,2)))) ;
    chi = any(chi,1)
    %}
    for j = 1:size(C,1)
        x0 = voxel_cor(1) - C(j,1);
        y0 = voxel_cor(2) - C(j,2);
        z0 = voxel_cor(3) - C(j,3);
        v = [cos(C(j,4)).*sin(C(j,5)) sin(C(j,4)).*sin(C(j,5)) cos(C(j,5))]';
        if (C(j,7) > acos( [x0 y0 z0]*v./(sqrt(sum(abs([x0 y0 z0]).^2,2)))));
            in(:,p) = M(j,:)';
            p = p + 1;
        end
    end

