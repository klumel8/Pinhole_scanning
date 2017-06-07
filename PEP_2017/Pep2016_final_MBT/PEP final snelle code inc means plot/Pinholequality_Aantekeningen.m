%% Propedeutic End Project      -   June 2014
% Hidde van Ulsen, Floris van der Gronden and Thijs van de Mortel
clear all; close all;
[x,y,z,theta,phi,alpha] = pinholes2;
% load Pinholes
%% System Specifications

Lx =20; Ly = 20; Lz = 20;               % Size of volume
dx = 2; dy = 2; dz = 2;           % Step size for voxels

gridrow = 19;                           % Gridsize: number of sections in the vertical direction
gridcol = 32;                           % Gridsize: number of sections in the horizontal direction

%% Voxel part

Nx = Lx/dx +1;
Ny = Ly/dy +1;
Nz = Lz/dz +1;
quality = zeros(Nx,Ny,Nz);  

Npin = length(x);
Ncircle = (gridrow-1)/2*gridcol +1;

rpin=zeros(Npin,3);

for l = 1:Npin
    rpin(l,:) = [x(l) y(l) z(l)];
end

normpin = zeros(Npin,3);
for n = 1:Npin
    normpin(n,:) = [-cos(phi(n))*sin(theta(n)) -sin(phi(n))*sin(theta(n)) cos(theta(n))];
end

xvox = (-Lx/2+(0:Nx-1)*dx);
yvox = (-Ly/2+(0:Ny-1)*dy);
zvox = (-Lz/2+(0:Nz-1)*dz);

%% Circles

% CircleGen input: gridrow/gridcol (= gridsize)
C = CircleGen(gridrow,gridcol);
% CircleGen output: 3D-matrix C with size (fx X fy X fz) with 
% fx = different coordinates on the circles, 
% fy = the actual x,y,z coordinates of the points on the circles
% fz = different circles


% Polarize input: 3D-matrix C
Cp = Polarize(C);
% Polarize output: 3D-matrix Cp which is the same as matrix C, but in
% stead of carthesian coordinates, you know have spherical coordinates. In
% the fy direction you now have 2 columns (first theta, then phi).

% CircleGrid input: Cp, gridrow, gridcol (= circles & gridsize)
gridcir=CircleGrid(Cp,gridrow,gridcol);       %(sinusvorm in grid)
% CirkelGrid output: 3D-matrix gridcir with in the first two dimensions the
% grid which 'displays' the circle line with 1's. In the third dimension
% are the different circles.



%% Determining Quality
for i = 1:Nx;
    for j = 1:Ny;
        for k = 1:Nz;
            FOV=zeros(Npin,1);
            Voxel=repmat([xvox(i) yvox(j) zvox(k)],[size(rpin,1),1]);
            rpinvox = (-Voxel + rpin);
            beta=acos((diag(-normpin*rpinvox'))./sqrt(sum(rpinvox.^2,2)));
            FOV(beta<=alpha/2)=1;
            FOV(beta>=pi-alpha/2)=1;
            %FOV geeft 1 als beta binnen de alfa (blauw) ligt
            
            phivox = atan2(rpinvox(:,2),(rpinvox(:,1)));
            thetavox = atan2(sqrt((rpinvox(:,1)).^2 + (rpinvox(:,2)).^2),-(rpinvox(:,3)));
                   
            phivox(phivox<0)=phivox(phivox<0)+2*pi;
            thetavox(thetavox<0)=thetavox(thetavox<0)+2*pi;
            
            
            grid = zeros(gridrow,gridcol);
            
            % MakeGridLoc input:  gridrow,gridcol,phi,theta,FOV           
            gridloc = MakeGridLoc(gridrow,gridcol,phivox,thetavox,FOV,Npin);
            % MaakGridloc output: 2D matrix gridloc. 
            % It is a 2D representation of the spherical surface seen as if standing inside the voxel.
            % It gives a 1 when there is a pinhole in that section and it
            % gives a 0 when there is no pinhole in that section.       
                        
                       
            tester=false(1,1,Ncircle);
            for m = 1:Npin
                if FOV(m,:) == 1
                    tester=tester | gridcir(gridloc(m,1),gridloc(m,2),:);
                end
            end
            
            quality(i,j,k) = sum(tester)/Ncircle;
        end
    end
end
temp=quality(1:11,1:11,1:11);  %pak een klein deel in het midden van het volume
meanplot = [mean(quality(:)) mean(temp(:)) median(quality(:)) min(temp(:)) median(temp(:)) max(quality(:))];