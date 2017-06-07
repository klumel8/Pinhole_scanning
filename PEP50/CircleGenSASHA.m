%function C = CircleGen(gridrow,gridcol)
%a = 0 %
Nc = 4*max(gridrow,gridcol);
dg = 2*pi/Nc;
ngrid = gridcol*gridrow;
C = zeros(Nc,3,ngrid); %
Res = zeros(Nc,2,ngrid);
tic
i = 0:gridcol-1;
phim = i*2*pi/gridcol;
j = 0:gridrow-1;
thetam = j*pi/2/gridrow;
a = 1:ngrid;
gamma = 1:Nc;
[Phi,Theta,Gamma] = meshgrid(phim,thetam,gamma);
Circles = cat(4,cat(4,times(times(cos(dg*Gamma),sin(Theta)),cos(Phi))-times(sin(dg*Gamma),sin(Phi)), times(times(cos(dg*Gamma),sin(Theta)),sin(Phi))+times(sin(dg*Gamma),cos(Phi))), times(cos(dg*Gamma),cos(Theta)));
C = reshape(Circles,[ngrid,Nc,3]);
C = permute(C,[2,3,1]);
Res(:,1,:) = acos(C(:,3,:));
Res(:,2,:) = atan2(C(:,2,:),C(:,1,:));
Res(:,2,:) = 2*pi*(Res(:,2,:)<0)+Res(:,2,:)-2*pi*(Res(:,2,:)>0);
toc
