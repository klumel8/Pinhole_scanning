%% Propedeutic End Project      -   June 2016
% Inge Faber, Wilbert Ras, Eva Spreeuw


close all; clc; clear all;
load VECTOR_collimator


%% Data unvolding 
x=VECTOR_collimator(:,1);               
y=VECTOR_collimator(:,2);               
z=VECTOR_collimator(:,3);          
theta=VECTOR_collimator(:,5);     
phi=VECTOR_collimator(:,4);
alpha=VECTOR_collimator(:,6);    
d=VECTOR_collimator(:,7); 

%% System Specifications

Lx = 25; Ly = 25; Lz = 80;              % Size of volume
dx = 1; dy = 1; dz = 1;                 % Step size for voxels

gridrow = 75;                           % Gridsize: number of sections in the vertical direction
gridcol = 150;                           % Gridsize: number of sections in the horizontal direction

radius_coll = 20; %mm                   %radius van colliminator
    %for ST
phantom_size = [Lx Ly Lz];
R_transaxial = 7;                       %raidus of spiral
N_Pos_rev = 4.5;                        %number of steps per revolution of spiral
N_Pos_tot = 25;                         %total number of steps

%% Voxel part

Nx = Lx/dx +1;
Ny = Ly/dy +1;
Nz = Lz/dz +1;
Uquality = zeros(Nx,Ny,Nz);  

Npin = length(x);
Ncircle = (gridrow-1)/2*gridcol +1;

rpin=zeros(Npin,3);

for l = 1:Npin
    rpin(l,:) = [x(l) y(l) z(l)];
end

normpin = zeros(Npin,3);
for n = 1:Npin
    normpin(n,:) = [cos(phi(n))*sin(theta(n)) sin(phi(n))*sin(theta(n)) cos(theta(n))];
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
gridcir=CircleGrid(Cp,gridrow,gridcol);
% CirkelGrid output: 3D-matrix gridcir with in the first two dimensions the
% grid which 'displays' the circle line with 1's. In the third dimension
% are the different circles.

gridloc = cell(1,Nz);
% relative movements due to bed positions
%Bed_Pos_helix input: phantom_size,R_transaxial,N_Pos_rev,N_Pos_tot (volume,radius, position per revolution, total number of positions).
movement = Bed_Pos_helix(phantom_size,R_transaxial,N_Pos_rev,N_Pos_tot);
%Bed_Pos_helix output: x, y and z arrays discribing ST.

%% Loop over total volume
for l=1:Nx;
     for m=1:Ny;
          for n=1:Nz;    
              gridloc{l}{m}{n} = [];
             
              for b = 1:N_Pos_tot;                          
                    FOV=zeros(Npin,1);
                    % Voxel position + relative movements due to bed positions
                    Voxel=repmat([xvox(l) yvox(m) zvox(n)]+movement(b,:),[size(rpin,1),1]);
                            
                     if  sqrt(Voxel(1,1)^2+Voxel(1,2)^2) <= radius_coll
                     rpinvox = (-Voxel + rpin);
                     
                     % Deterimining if voxel is in FOV
                     beta=acos(sum(-normpin.*rpinvox,2)./sqrt(sum(rpinvox.^2,2)));
                     FOV(beta<=alpha/2)=1;
                     FOV(beta>=pi-alpha/2)=1;
 
                     phivox = atan2(rpinvox(:,2),(rpinvox(:,1)));
                     thetavox = atan2(sqrt((rpinvox(:,1)).^2 + (rpinvox(:,2)).^2),-(rpinvox(:,3)));
 
                     phivox(phivox<0)=phivox(phivox<0)+2*pi;
                     thetavox(thetavox<0)=thetavox(thetavox<0)+2*pi;

                     % MakeGridLoc input:  gridrow,gridcol,phi,theta,FOV           
                     gridloc{l}{m}{n} = cat(1,gridloc{l}{m}{n},MakeGridLoc(gridrow,gridcol,phivox,thetavox,FOV));
                     % MaakGridloc output: 2D matrix gridloc. 
                     % It is a 2D representation of the spherical surface seen as if standing inside the voxel.
                     % It gives a 1 when there is a pinhole in that section and it
                     % gives a 0 when there is no pinhole in that section.       
                            
                     end
              end
              %%  gridsoptellen
              TOTgridloc = gridloc{l}{m}{n};
              %%  quality uitrekenen
              tester = false(1,1,Ncircle);
              sizeTOTgridloc = size(TOTgridloc);
              for o = 1:sizeTOTgridloc(1);
                  if TOTgridloc(o,3) == 1
                  tester=tester | gridcir(TOTgridloc(o,1),TOTgridloc(o,2),:);
                  end
              end

              Vquality(l,m,n) = sum(tester)/Ncircle;

          end
      end
end






