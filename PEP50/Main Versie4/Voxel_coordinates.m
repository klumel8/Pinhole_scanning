
function [ voxels ] = Voxel_coordinates( Lx, Ly, Lz, dx, dy ,dz )

%voxel coordinaten Geeft de coordinaten van de voxels  terug 

% the rows is de y index van 
% the colums is de x index 
% the 3d component is the z index

Nx = Lx/dx;
Ny = Ly/dy;
Nz = Lz/dz;

voxels = zeros(Nx*Ny*Nz,3);
xvox = (-Lx/2+(0:Nx-1)*dx);
yvox = (-Ly/2+(0:Ny-1)*dy);
zvox = (-Lz/2+(0:Nz-1)*dz);

i = 0;

for l=1:Nx;
     for m=1:Ny;
      for k = 1:Nz;
          i = i+ 1;
       voxels(i,:) = [xvox(l),yvox(m),zvox(k)];      
      end    
     end
end

end