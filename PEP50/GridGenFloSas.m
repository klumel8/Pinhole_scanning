clear all

gridrow = 129; %Number of rows
gridcol = 256; %Number of colms
halvedrow = (gridrow-1)/2+1;

rowleng=pi/halvedrow/2; %Height of a row
colleng=2.*pi/gridcol; %Width of a col

Nc = max(halvedrow,gridcol)*5; %Number of steps of gamma
dg = 2*pi/Nc;  %Stepssize of gamma parameter
ngrid = gridcol*halvedrow;% gridcol*gridrow; %Number of circles  
gridcir=zeros(ngrid,gridrow,gridcol); %Set the ngrid circlegrids on false

i = 0:gridcol-1; 
phim = i*2*pi/gridcol; %Vector from 0 to 2pi-colleng in  ; so phi is a the row (x)

j = 0:halvedrow-1; %j = 0:gridrow-1; 
thetam = j*pi/2/gridrow; %Vector from 0 to pi- rowleng ; So theta is the collumn (y)

a = 1:ngrid; 
gamma = 1:Nc; 
[Phi,Theta,Gamma] = meshgrid(phim,thetam,gamma);

Circles = cat(4,cat(4,times(times(cos(dg*Gamma),sin(Theta)),cos(Phi))-times(sin(dg*Gamma),sin(Phi)), times(times(cos(dg*Gamma),sin(Theta)),sin(Phi))+times(sin(dg*Gamma),cos(Phi))), times(cos(dg*Gamma),cos(Theta)));
C = reshape(Circles,[ngrid,Nc,3]);
C = permute(C,[2,3,1]);
Res(:,1,:) = acos(C(:,3,:)); %Theta 
Res(:,2,:) = atan2(C(:,2,:),C(:,1,:)); %Phi
Res(:,2,:) = 2*pi*(Res(:,2,:)<0)+Res(:,2,:)-2*pi*(Res(:,2,:)>2*pi);

index(:,1,:) = floor(Res(:,1,:)/rowleng)+1;%Theta
index(:,2,:) = floor(Res(:,2,:)/colleng)+1;%Phi

linearIndex = zeros(Nc,1,ngrid);
linearIndex = halvedrow*index(:,2,:)+index(:,1,:);

gridcir(a,linearIndex(:))
%{
weight=zeros(ngrid,1);
tic
grid = zeros(gridrow,gridcol);
for i = 1:length(index(:,:,1))
% NOTE NOTE NOTE zou dit niet |grid(index(i,1,j),index(i,2,j)) =1;| moeten zijn waarin we ook nog eens door j loopen voor elk punt(1 tot Nc)?
    grid(index(i,1,1),index(i,2,1)) =1;
end
toc
%}
