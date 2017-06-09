clear all

gridrow = 201; %Number of rows
gridcol = 400; %Number of colms

rowleng=pi/gridrow; %Height of a row
colleng=2*pi/gridcol; %Width of a col

numberSteps = max(gridrow,gridcol)*5; %Number of steps of gamma
dg = 2*pi/numberSteps;  %Stepssize of gamma parameter
numberCircles = gridcol*((gridrow-1)/2+1);% gridcol*gridrow; %Number of circles  
gridcir=false(numberCircles,gridrow,gridcol); %Set the ngrid circlegrids on false

i = 0:gridcol-1; 
phim = i*2*pi/gridcol; %Vector from 0 to 2pi-colleng in  ; so phi is a the row (x)

j = 0:(gridrow-1)/2; %j = 0:gridrow-1; 
thetam = j*pi/2/gridrow; %Vector from 0 to pi- rowleng ; So theta is the collumn (y)

a = 1:numberCircles; 
gamma = 1:numberSteps; 
[Phi,Theta,Gamma] = meshgrid(phim,thetam,gamma);

Circles = cat(4,cat(4,...
    times(times(cos(dg*Gamma),sin(Theta)),cos(Phi))-times(sin(dg*Gamma),sin(Phi)),...
    times(times(cos(dg*Gamma),sin(Theta)),sin(Phi))+times(sin(dg*Gamma),cos(Phi))), ...
    times(cos(dg*Gamma),cos(Theta)));
C = reshape(Circles,[numberCircles,numberSteps,3]);
C = permute(C,[2,3,1]);
Res(:,1,:) = acos(C(:,3,:)); %Theta 
Res(:,2,:) = atan2(C(:,2,:),C(:,1,:)); %Phi
Res(:,2,:) = 2*pi*(Res(:,2,:)<0)+Res(:,2,:)-2*pi*(Res(:,2,:)>2*pi);
%tic
index(:,1,:) = floor(Res(:,1,:)/rowleng);%Theta
index(:,2,:) = floor(Res(:,2,:)/colleng);%Phi


linearIndex = gridrow*index(:,2,:)+index(:,1,:);
linearIndex = squeeze(linearIndex);

labels= numberCircles*meshgrid(a-1,gamma);
linearIndex = linearIndex+labels;
linearIndex = reshape(linearIndex,[numberSteps*numberCircles,1]);



gridcir(linearIndex(:,1))=1;
%toc

%weight=zeros(ngrid,1);
%{
tic
grid = zeros(gridrow,gridcol);
for i = 1:length(index(:,:,1))
    grid(index(i,1,1),index(i,2,1)) =1;
end
toc
%}
