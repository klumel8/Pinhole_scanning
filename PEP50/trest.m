M = [x y z phi -theta d alpha];

quiver3(M(:,1),M(:,2),M(:,3),cos(M(:,4)).*sin(M(:,5)),sin(M(:,4)).*sin(M(:,5)),cos(M(:,5)));
hold on 
plot3(x_voxel,y_voxel,z_voxel,'.')

