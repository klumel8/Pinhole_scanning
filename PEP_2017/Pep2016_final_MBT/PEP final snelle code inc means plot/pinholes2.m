function [X,Y,Z,theta,phi,alpha]=pinholes2()

z=[-50;-30;0;30;50];
x=[-50;-30;0;30;50];
ztheta=[50;73;90;107;130]/180*pi+pi;
xtheta=[50;73;90;107;130]/180*pi+pi;
n_x=length(x);
n_z=length(z);
[Z,X]=meshgrid(x,z);
Y=30*ones(size(X));
X=[-X(:);-X(:)];
Y=[Y(:);-Y(:)];
Z=[Z(:);Z(:)];
[kaderaxis_yz,kaderaxis_xy]=meshgrid(ztheta,xtheta);

%figure(10);clf;
%drawslits(n_z,z,2,15/180*pi*ones(size(z)),ztheta,5,110,0,55, 400,120,5)
%figure(11);clf;
%drawslits(n_x,x,2,15/180*pi*ones(size(z)),xtheta,5,110,0,55, 500,120,5)

axisph=[1./tan(kaderaxis_xy(:)');ones(size(kaderaxis_xy(:)'));1./tan(kaderaxis_yz(:)')];
tmp=[1./tan(kaderaxis_xy(:)');-ones(size(kaderaxis_xy(:)'));1./tan(kaderaxis_yz(:)')];
axisph=[axisph,tmp];clear tmp;%add second set of pinholes
for i=1:50
    axisph(:,i)=axisph(:,i)/norm(axisph(:,i));
end;

theta=acos(axisph(3,:));%- for flipping z-axis
phi=atan2(axisph(2,:),axisph(1,:));
for i = 1:50
    if phi(i) < 0
        phi(i) = phi(i) + 2* pi;
    end
end
alpha = 0.5632.*ones(50,1);

%normal=@(theta,phi) [cos(phi).*sin(theta);sin(phi).*sin(theta);cos(theta)];
%N=normal(theta,phi);
%figure(20);quiver3(X,Y,Z,N(1,:)',N(2,:)',N(3,:)');