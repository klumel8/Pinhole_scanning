%% Propedeutic End Project      -   June 2016
% Rik, Joram, Willemijn: bewerking van: 
% Hidde van Ulsen, Floris van der Gronden and Thijs van de Mortel
clc; clear all;
%for-loop voor verschillende afstanden

L=40;
s=6;
L0=L-s;


width=150; % x-richting
chesttonipple=110; %z-richting
thickness=55; %y-richting

nalpha = 1;
amin = 40;
amax = 40;
alpha_NEW = linspace(amin,amax,nalpha);
afstanden = 2.*L0.*tan(alpha_NEW.*pi/180); %maak afstanden zo dat alpha_yz=alpha_NEW (design)

time=tic;
meanplot=zeros(nalpha,7);

R_target=6.0;
do_seq=1;
Ri=3.2;
d_factor=1.0;

%variabelen voor BET-posities
xsmin=-75;
xsmax=75;
sgx=6;                              %(xsmax-xsmin)/nxstap; %sgx= stapgrootte x
nxstap=ceil((xsmax-xsmin)/sgx);     %nxstap= aantal stappen x
sgz=2;                              %sgz = stapgrootte z
%% System Specifications

Lx =150; Ly = 54; Lz = 110;            % Size of volume
dx =2; dy =2; dz =2;                   % Step size for voxels

gridrow = 129;                          % Gridsize: number of sections in the vertical direction
gridcol = 256;                        % Gridsize: number of sections in the horizontal direction

%% Voxel part

Nx = Lx/dx +1;
Ny = Ly/dy +1;
Nz = Lz/dz +1;
Ncircle = (gridrow-1)/2*gridcol +1;
xvox = (-Lx/2+(0:Nx-1)*dx); 
yvox = (-Ly/2+(0:Ny-1)*dy);
zvox = (0:Lz);
%zorgt ervoor dat het midden van het volume in de origin ligt. xvox geeft
%alle plaatsen van de pinholes in de x-richting
%% Circles

% C = CircleGen(gridrow,gridcol);
% % CircleGen output: 3D-matrix C with size (fx X fy X fz) with
% % fx = different coordinates on the circles,
% % fy = the actual x,y,z coordinates of the points on the circles
% % fz = different circles
% % Polarize input: 3D-matrix C
% Cp = Polarize(C);
% % Polarize output: 3D-matrix Cp which is the same as matrix C, but in
% % stead of carthesian coordinates, you know have spherical coordinates. In
% % the fy direction you now have 2 columns (first theta, then phi).
% % CircleGrid input: Cp, gridrow, gridcol (= circles & gridsize)
% gridcir=CircleGrid(Cp,gridrow,gridcol);

[gridcir2,~]=CircleGrid2(gridrow,gridcol);


for ii=1:length(afstanden)
    %fprintf('Afstand: %d/%d\n',ii,length(afstanden));
    afstand=afstanden(ii); 

[Sens,d,nph,S,rph,axisph,alpha_yz,nzstap]=design_201307_short_det(xsmin,xsmax,nxstap,sgz,afstand,L,R_target,do_seq,Ri,d_factor);

x=rph(1,:);
y=rph(2,:); 
w=rph(3,:);

if  do_seq>0
        for ji=1:(nxstap-1) 
            xa=rph(1,:)+xsmin+ji*sgx*ones(1,length(rph)); 
            x=[x xa];
            y=[y rph(2,:)];
            w=[w rph(3,:)];
        end
    xb=x;
    yb=y;
    wb=w;
    for wp=1:(nzstap-1)
             wa=wb-wp*sgz*ones(1,length(wb));
             w=[w wa];
             x=[x xb];
             y=[y yb];
    end
end
z=685-w;
rph=[x; y; z];
alpha= alpha_yz;

figure(1); %hold on
plot3(rph(1,:),rph(2,:),rph(3,:),'.')
xlabel('x-as'); ylabel('y-as'); zlabel('z-as');
view(3);

quality = zeros(Nx,Ny,Nz);  
Npin = length(x);

rpin=zeros(Npin,3);

for l = 1:Npin
    rpin(l,:) = [x(l) y(l) z(l)];
end


normpin=axisph; % aanpassing voor design-201307
Lnp=length(normpin);
normpin2=zeros(3,2*Lnp);
for irik=1:Lnp;
    normpin2(:,irik)=([-1 0 0; 0 -1 0; 0 0 1]*axisph(:,irik));
    normpin2(:,irik+Lnp)=([-1 0 0; 0 1 0; 0 0 1]*axisph(:,irik));
end

normpin3=normpin2;
for i=1:((nxstap*nzstap)-1)
    normpin3=[normpin3 normpin2];
end

% normpin = zeros(Npin,3);
% for n = 1:Npin
%     normpin(n,:) = [-cos(phi(n))*sin(theta(n)) -sin(phi(n))*sin(theta(n)) cos(theta(n))];
% end




breast=false(Nx,Ny,Nz);

%fprintf('Setup Done. Starting Quality:\n');
%% Determining Quality
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            FOV=false(Npin,1);
            Voxel=repmat([xvox(i) yvox(j) zvox(k)],[size(rpin,1),1]);
            rpinvox = (-Voxel + rpin);
            beta=acos((sum(-normpin3'.*rpinvox,2))./sqrt(sum(rpinvox.^2,2)));
            FOV(beta<=alpha/2)=true;
            FOV(beta>=pi-alpha/2)=true;
                %dit geeft een 1 als de voxel in the field of view van de
                %pinhole zit en 0 als de voxel er niet in zit
                        
            phivox = atan2(rpinvox(:,2),(rpinvox(:,1)));
            thetavox = atan2(sqrt((rpinvox(:,1)).^2 + (rpinvox(:,2)).^2),-(rpinvox(:,3)));
                   
            phivox(phivox<0)=phivox(phivox<0)+2*pi;
            thetavox(thetavox<0)=thetavox(thetavox<0)+2*pi;
            
            phivox=phivox(FOV);
            thetavox=thetavox(FOV);
            
            grid = zeros(gridrow,gridcol);
            
            % MakeGridLoc input:  gridrow,gridcol,phi,theta,FOV           
            gridloc = MakeGridLoc(gridrow,gridcol,phivox,thetavox);
            % MaakGridloc output: 2D matrix gridloc. 
            % It is a 2D representation of the spherical surface seen as if standing inside the voxel.
            % It gives a 1 when there is a pinhole in that section and it
            % gives a 0 when there is no pinhole in that section.       
                              
            tester=false(Ncircle,1);
            
            for m = 1:length(phivox)
                tester=tester | gridcir2(:,gridloc(m,1),gridloc(m,2));
            end
            
            quality(i,j,k) = sum(tester)/Ncircle;
                          if((xvox(i)/(width/2))^2+((zvox(k))/chesttonipple)^2)<=1 && zvox(k)>10 && yvox(j)>-thickness/2&& yvox(j)<thickness/2
                              breast(i,j,k)=true;
                         end
        end
    end
end
 qualitymask=quality(breast);
    %qualitymask=0;

temp=quality;%temp=quality(1:11,1:11,1:11); 
meanplot(ii,:) = [mean(quality(:)) mean(temp(:)) median(quality(:)) min(temp(:)) median(temp(:)) max(quality(:)) mean(qualitymask(:))];
%{
imwrite(squeeze(quality(45,:,:))',sprintf('fig1_%02d.png',ii));
imwrite(squeeze(quality(:,15,:))',sprintf('fig2_%02d.png',ii));
imwrite(squeeze(quality(:,:,62)'),sprintf('fig3_%02d.png',ii));
imwrite(squeeze(quality(:,:,72)'),sprintf('fig4_%02d.png',ii));
%}


%save(fn);
end
toc(time)
figure(2); hold on
plot(alpha_NEW,meanplot(:,1),'r');
plot(alpha_NEW,meanplot(:,7),'b');
%save('quality.mat');