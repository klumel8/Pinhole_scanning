function [Sens,d,nph,S,rph,axisph,alpha_yz,nzstap]=design_201307_short_det(xsmin,xsmax,nxstap,sgz,afstand,L,R_target,do_seq,Ri,d_factor)
% calculates sensitivity and resolution in a subvolume using
% offsets for pinholes.
% L = distance detector (typical 40 to 150 mm)(afstand tussen middel van volume en detector) 
% R_target = designed resolution, say 6.0 mm
% do_seq says if bed positions need to be modelled, 0=off,1=whole breast scan, 2 = focussed scan
% Ri = Detector resolution, leave at 3.2 mm
% d_factor is for internal use, do not set, or set to 1.0
tic;
if nargin<5
    d_factor=1.0;
end;
%%%%%%%%%START DESIGN COLLIMATOR%%%%%%%%%%%

s=6; %distance pinhole from breast (afstand tussen ph en het midden van volume)
dist=0';
L0=L-s;  %(L0 is afstand tussen detector en ph)
%% side view
DX=140*1.072;%Detector size
DY=234*1.072;
% FZ=(s*DX + L0*maxchestnipple + DX*maxthickness)/(DX - maxchestnipple);
% if (isinf(FZ))
%     FZ=1e16;
% end;
FZ=1e14; 
fp=[-(FZ-27.5),685];% focal distance

pivot=[27.5+s+L0,685-DX]; %focal distance
alpha=0;

%detector intercept functions

detectoryz=@(y,z) ([y,z]-fp)*(dot(pivot-fp,[-cos(alpha);sin(alpha)])/dot([y,z]-fp,[-cos(alpha);sin(alpha)]))+fp;
directionyz=@(y,z) ([y,z]-fp)/norm([y,z]-fp);
detectoryz2=@(dir,point) dir*(dot(pivot-point,[-cos(alpha);sin(alpha)])/dot(dir,[-cos(alpha);sin(alpha)]))+point;

z_front=685;
y_front=27.5+s+dist;
y_target=27.5;

%define some other constants
mu_=3.3892;

dist=s;
%determining pinhole diameter
d=sqrt(((L0)^2*R_target^2-(y_target+dist)^2*Ri^2)/((y_target+dist)^2+2*(y_target+dist)*(L0)+(L0)^2))*d_factor;
if d<0.1 || imag(d)>0
    Sens=[0 0 0 0];
    d=0;
    nph=0;
    return;
end;
dist=0;

%rotation matrices
alpha_yz=atan((afstand)/2/L0);
alpha_xy=atan((afstand)/2/L0);
R_yz_up=[cos(alpha_yz),-sin(alpha_yz);sin(alpha_yz),cos(alpha_yz)];
R_yz_down=[cos(-alpha_yz),-sin(-alpha_yz);sin(-alpha_yz),cos(-alpha_yz)];
R_xy_left=[cos(alpha_xy),-sin(alpha_xy);sin(alpha_xy),cos(alpha_xy)];
R_xy_right=[cos(-alpha_xy),-sin(-alpha_xy);sin(-alpha_xy),cos(-alpha_xy)];




rdetector=detectoryz(y_front,z_front);

clear rph axiskader X Y Z;
r_lower=[rdetector(1),z_front];
nphz=0;
while ((r_lower(2)>pivot(2))) && (nphz<100)
    %find location of next pinhole
    r_upper=[inf,inf];i=0;
    while (r_upper(2)>r_lower(2)-3) && (i<5000) % deze while-loop zorgt ervoor dat voor alle z-coordinaten pinholes worden gemaakt
        z_front=z_front-0.05;
        dir=directionyz(y_front,z_front);
        upper_dir=R_yz_up*(dir');
        r_upper=detectoryz2(upper_dir',[y_front,z_front+d/2]);
        i=i+1;
    end;
    %draw upper septawall line
    Z1(nphz+1)=z_front;
    rdetector=detectoryz(y_front,z_front);
    kaderaxis_yz1(nphz+1)=atan2(rdetector(2)-fp(2),rdetector(1)-fp(1));
    kader_opening_yz_upper1(nphz+1)=alpha_yz;
    kader_opening_yz_lower1(nphz+1)=alpha_yz;
    
    %determine location of lowest detector image
    dir=directionyz(y_front,z_front);
    lower_dir=R_yz_down*(dir');
    r_lower=detectoryz2(lower_dir',[y_front,z_front-d/2]);
    det_lower1(nphz+1)=r_lower(2);
    det_upper1(nphz+1)=r_upper(2);
    nphz=nphz+1;
    
    if (lower_dir(1)<sin(10/180*pi))
        fprintf('Breaking z placement due to looking backward\n');
        break;
    end;
    if r_lower(2)>r_upper(2)
        fprintf('Upper edge detected below of lower edge. Stop');
        Sens=[0 0];
        return;
    end;
    % fprintf('Upper: %g\t\tLower: %g\n',r_upper(2),r_lower(2));
    
end
%det=[0,-1];
nph=0;



for iz=1:nphz %deze for-loop zorgt ervoor dat voor ALLE z-coordinaten van pinholes naar links en naar rechts toe (in de x-richting) de bijbehorende pinholes worden gemaakt.--> hele raster aan pinholes
    %get important statistics from z axis
    z_front=Z1(iz);
    rdetector=detectoryz(y_front,z_front);
    angle=kaderaxis_yz1(iz);
    Det_dist(iz)=pivot(1);
    
    fp=[0,-40+27.5+s];
    detector=@(x,y) ([x,y]-fp)*(dot([0,Det_dist(iz)]-fp,[0;-1])/dot([x,y]-fp,[0;-1]))+fp;
    directionxy=@(x,y) ([x,y]-fp)/norm([x,y]-fp);
    detectorxy2=@(dir,point) dir*(dot([0,Det_dist(iz)]-point,[0;-1])/dot(dir,[0;-1]))+point;
    
    x_front=0;
    y_front=27.5+s+dist;
    %add cones
    %init vars
    rdetector=detector(x_front,y_front);
    
    %Save middle pinhole
    X(nph+1)=0;
    Z(nph+1)=Z1(iz);
    kaderaxis_xy(nph+1)=atan2((rdetector(2)-y_front),rdetector(1)-x_front);
    kaderaxis_yz(nph+1)=kaderaxis_yz1(iz);
    kader_opening_yz_lower(nph+1)=kader_opening_yz_lower1(iz);
    kader_opening_yz_upper(nph+1)=kader_opening_yz_upper1(iz);
    kader_opening_xy_l(nph+1)=alpha_xy;
    kader_opening_xy_r(nph+1)=alpha_xy;
    det_lower(nph+1)=det_lower1(iz);
    det_upper(nph+1)=det_upper1(iz);
    
    %determine location of rightest detector image
    dir=directionxy(x_front,y_front);
    right_dir=R_xy_right*(dir');
    r_right=detectorxy2(right_dir',[x_front+d/2,y_front]);
    det_right(nph+1)=r_right(1);
    left_dir=R_xy_left*(dir');
    r_left=detectorxy2(left_dir',[x_front-d/2,y_front]);
    det_left(nph+1)=r_left(1);
    clear left_dir r_left;
    
    %move coordinates to next ph 
    nph=nph+1;
    k=0;
    %right side 
    while ((r_right(1)<DY/2+0)) && (k<50)
        %find location of next pinhole
        r_left=[0,0];i=0;
        while (r_left(1)<r_right(1)+3) && (i<5000)
            x_front=x_front+0.02;
            dir=directionxy(x_front,y_front);
            left_dir=R_xy_left*(dir');
            dir2=[dir(2);-dir(1)];
            
            r_left=detectorxy2(left_dir',[x_front-d/2*dir2(1),y_front-d/2*dir2(2)]);
            i=i+1;
        end;
        if (r_left(1)<r_right(1)+3)
            printf('Stopped placing pinholes due to exeeding max #iter.');
        end;
        X(nph+1)=x_front;
        rdetector=detector(x_front,y_front);
        kaderaxis_xy(nph+1)=atan2((rdetector(2)-y_front),rdetector(1)-x_front);
        kader_opening_xy_r(nph+1)=alpha_xy;
        kader_opening_xy_l(nph+1)=alpha_xy;
        
        %Save z values
        Z(nph+1)=Z1(iz);
        kaderaxis_yz(nph+1)=kaderaxis_yz1(iz);
        kader_opening_yz_lower(nph+1)=kader_opening_yz_lower1(iz);
        kader_opening_yz_upper(nph+1)=kader_opening_yz_upper1(iz);
        %determine location of lowest detector image
        det_lower(nph+1)=det_lower1(iz);
        det_upper(nph+1)=det_upper1(iz);
        
        dir=directionxy(x_front,y_front);
        right_dir=R_xy_right*(dir');
        dir2=[dir(2);-dir(1)];
        
        r_right=detectorxy2(right_dir',[x_front+d/2*dir2(1),y_front+d/2*dir2(2)]);
        det_right(nph+1)=r_right(1);
        det_left(nph+1)=r_left(1);
        
        
        nph=nph+1;
        
        if right_dir(2)<sin(10/180*pi)
            fprintf('Breaking x right placement due to looking backward\n');
            break; %pinhole looks partly upward.
        end;
        if r_right(1)<r_left(1)
            fprintf('Right edge detected left of left edge. Stop');
            Sens=[0 0 0];
            return;
        end;
        %fprintf('Left:\t\t%g\t\tright:\t%g\n',r_left(1),r_right(1));
        
        k=k+1;
    end
    
    
    %left side
    x_front=0;
    %determine location of rightest detector image
    dir=directionxy(x_front,y_front);
    left_dir=R_xy_left*(dir');
    r_left=detectorxy2(left_dir',[x_front-d/2,y_front]);
    k=0;
    while ((r_left(1)>-DY/2-0)) && (k<50)
        %find location of next pinhole
        r_right=[inf,inf];i=0;
        while (r_right(1)>r_left(1)-3) && (i<5000)
            x_front=x_front-0.02;
            dir=directionxy(x_front,y_front);
            dir2=[dir(2);-dir(1)];
            right_dir=R_xy_right*(dir');
            r_right=detectorxy2(right_dir',[x_front+d/2*dir2(1),y_front+d/2*dir2(2)]);
            i=i+1;
        end;
        if (r_right(1)>r_left(1)-3)
            printf('Stopped placing pinholes due to exeeding max #iter.');
        end;
        
        X(nph+1)=x_front;
        rdetector=detector(x_front,y_front);
        kaderaxis_xy(nph+1)=atan2((rdetector(2)-y_front),rdetector(1)-x_front);
        kader_opening_xy_r(nph+1)=alpha_xy;
        kader_opening_xy_l(nph+1)=alpha_xy;
        
        %Save z values
        Z(nph+1)=Z1(iz);
        kaderaxis_yz(nph+1)=kaderaxis_yz1(iz);
        kader_opening_yz_lower(nph+1)=kader_opening_yz_lower1(iz);
        kader_opening_yz_upper(nph+1)=kader_opening_yz_upper1(iz);
        %determine location of lowest detector image
        
        det_lower(nph+1)=det_lower1(iz);
        det_upper(nph+1)=det_upper1(iz);
        
        dir=directionxy(x_front,y_front);
        dir2=[dir(2);-dir(1)];
        
        left_dir=R_xy_left*(dir');
        r_left=detectorxy2(left_dir',[x_front-d/2*dir2(1),y_front-d/2*dir2(2)]);
        det_right(nph+1)=r_right(1);
        det_left(nph+1)=r_left(1);
        
        nph=nph+1;
        %fprintf('Left:\t\t%g\t\tright:\t%g\n',r_left(1),r_right(1));
        if left_dir(2)<sin(10/180*pi)
            fprintf('Breaking x left placement due to looking backward\n');
            break; %pinhole looks partly upward.
        end;
        if r_right(1)<r_left(1)
            fprintf('Right edge detected left of left edge. Stop');
            Sens=[0 0 0];
            return;
        end;
        
        
    end
    k=k+1;
end;

%%PLOTTING PINHOLES FOR 1 COLLIMATOR
%sideview
%figure(10);
% clf;
% drawslits(nphz,Z1-685+5,d,kader_opening_yz_lower1,pi/2+kaderaxis_yz1,s,L0-dist,dist,55,DX,-110,0);
%topview
%figure(11)
% clf;
% drawslits(nph/nphz,X,d,kader_opening_xy_r,-kaderaxis_xy,s,L0-dist,dist,55,DY,150,0);

%return

%Define pinholes

Y=ones(size(X))*(27.5+s+dist);% Y coordinates calulated from intersect of lowest and highest pinhole views.

%create list of location vectors of pinholes
rph=[X(:)';Y(:)';Z(:)'];
rph2=[X(:)';-Y(:)';Z(:)']; %zelf toegevoegd
%create axiskader
axiskader=[tan(pi/2-kaderaxis_xy(:)');ones(size(kaderaxis_xy(:)'));1./tan(pi/2-kaderaxis_yz(:)')];

% add pinholes on second detector second detector
rph=[rph,[rph(1,:);-rph(2,:);rph(3,:)]];
axiskader=[axiskader,[axiskader(1,:);-axiskader(2,:);axiskader(3,:)]];
tmp=kader_opening_yz_upper(:)';
kader_opening_yz_upper=[tmp,kader_opening_yz_lower(:)'];
kader_opening_yz_lower=[kader_opening_yz_lower(:)',tmp];
kader_opening_xy_left=[kader_opening_xy_l(:)',kader_opening_xy_r(:)'];
kader_opening_xy_right=[kader_opening_xy_r(:)',kader_opening_xy_l(:)'];
tmp=det_upper;
det_upper=[tmp,det_lower];
det_lower=[det_lower,tmp];
tmp=det_right;
det_right=[tmp,det_left];
det_left=[det_left,tmp];
% tot dit moment wel pinholes van beide platen meegenomen.
%number of pinholes
nph=length(rph)/2; % gedeeld door 2 omdat dit alleen pinholes zijn van ene plaat. 
%create axis for pinholes (unit vectors for locations);
% normalize axiskader vector
axisph=[1./tan(kaderaxis_xy(:)');ones(size(kaderaxis_xy(:)'));1./tan(pi/2-kaderaxis_yz(:)')];

for i=1:nph
    axiskader(:,i)=axiskader(:,i)/norm(axiskader(:,i));
    axisph(:,i)=axisph(:,i)/norm(axisph(:,i));
end;
%create pinhole openings
alpha_angle=ones(1,nph)*(2*atan(sqrt(alpha_xy^2 + alpha_yz^2)));
alpha=0;
%detector specification
det_width=DY;%width of detector active area;


%fprintf('done\nStarting calculation...\n');%end of config statement
%%DEBUG


%%%END DESIGN COLLIMATOR %%%%%%%


%%%%BEGIN SENSITIVITY AND RESOLUTION CALCULATION

nph=length(rph)/2;
r0det=[0;pivot(1);pivot(2)];
ndet=[0;-cos(alpha);sin(alpha)];

%Define volume for calculating resolution / sensitivity
xmin=-80;xmax=80;
ymin=-27.5;ymax=27.5;
zmin=685-125;zmax=685;
%zmin=692.8-1.5*50;zmax=692.8;
%xmin=-50;xmax=50;

%stepsizes -> grootte van voxels
dx=4;
dy=4;
dz=4;

%number of points along each dimension
nx=floor((xmax-xmin)/dx)+1;
ny=floor((ymax-ymin)/dy)+1;
nz=floor((zmax-zmin)/dz)+1;

%create position vectors
x=xmin+dx*(0:nx-1);
y=ymin+dy*(0:ny-1);
z=zmin+dz*(0:nz-1);

%Create meshgrids
[Y,X,Z]=meshgrid(y,x,z);

%allocate output arrays
S=zeros(nx,ny,nz);
%Covered=zeros(nx,ny,nz);
Rt=zeros(nx,ny,nz);

%keep track of intercept with detector
%Det_hitx=zeros(nph,nx,ny,nz);
%Det_hity=zeros(nph,nx,ny,nz);
%Det_hitz=zeros(nph,nx,ny,nz);
%dhx=[];dhy=[];dhz=[];


%calculate distance z
unique_z=unique(rph(3,:));
if numel(unique_z)==1
    distance_z=unique_z(1);
else
    distance_z=unique_z(2)-unique_z(1);
end

nzstap=ceil(distance_z/sgz);   %nzstap = aanstal stappen z-richting

seq_z=[];
if (do_seq==1)
    for i=1:nzstap %i=1:12;
            seq_z(i)=sgz*(i-1);
    end
    %seq_x=-78:8:78;
    seq_x=linspace(xsmin,xsmax,nxstap);
elseif (do_seq==2)
    for i=1:12;
       seq_z(i)=2*(i-1);
    end
    seq_x=[-10,-6,-4,-2,2,4,6,10];
else
    seq_x=0;
    seq_z=0;
end;

n_seq_x=length(seq_x);
%
n_seq_z=length(seq_z);

DeltaLk=log(2)/mu_;


% for iz=1:nz
%     
%     %fprintf('\b\b\b\b\b\b\b%3i/%3i',iy,ny);
%     %for ix=1:nx
%     for iy=1:ny
%         
%         for ix=1:nx
%             rsource=[x(ix);y(iy);z(iz)];
%             precalc=ndet/((r0det-rsource)'*ndet);
%             
%             for iseq=1:n_seq_x
%                 xseq=seq_x(iseq);
%                 
%                 for iph=1:nph
%                     rpinh=rph(:,iph);
%                     rpinh(1)=rpinh(1)+xseq;
%                     %find coordinates of this point
%                     
%                     
%                     v=rpinh-rsource;
%                     rdet=v/(v'*precalc) + rsource;
%                     if rdet(3) < det_lower(iph) || rdet(3) > det_upper(iph)  || rdet(1)<det_left(iph)+xseq || rdet(1)>det_right(iph)+xseq ||  ~(abs(rdet(1)-xseq)<det_width/2) || rdet(3)<pivot(2) || rdet(3) >685
%                         continue;% outside of active detector area, skipping
%                     end;
%                     
%                     dr=rsource-rpinh;
%                     
%                     h=abs(axisph(:,iph)'*dr);
%                     theta=asin(h/norm(dr));
%                     alpha=alpha_angle(iph);
%                     % the sensitivity function from Metzler et al, IEEE TMI 20(8) 2001, eq. (37)
%                     s=d^2*sin(theta)^3/(16*h^2)+sin(theta)^5*tan(alpha/2)^2/(8*h^2*mu_^2)*sqrt(1-cot(theta)^2/tan(alpha/2)^2)*...
%                         (1-cot(theta)^2/(tan(alpha/2)^2)+mu_*d*csc(theta)*cot(alpha/2));
%                     s=real(s);
%                     S(ix,iy,iz)=S(ix,iy,iz)+s;
%                     %Covered(ix,iy,iz)=Covered(ix,iy,iz)+1;
%                     
%                     
%                     dr2=rpinh-rdet;
%                     M=norm(dr2)/norm(dr);
%                     
%                     % effective diametr Metzler & Accorsi PMB 50 (2005) 5005-5017
%                     
%                     deff_par=d+DeltaLk*(tan(alpha/2)^2-cot(theta)^2)*cot(alpha/2)*sin(theta);
%                     deff_perp=sqrt((d+DeltaLk*tan(alpha/2)*sin(theta))^2-DeltaLk^2*cos(theta)^2);
%                     g_par=deff_par*(1+1/M);
%                     g_perp=deff_perp*(1+1/M);
%                     Rt_par=sqrt(g_par^2+Ri^2/M^2);
%                     Rt_perp=sqrt(g_perp^2+Ri^2/M^2);
%                     R=(Rt_par+Rt_perp)/2;
%                     R=real(R);
%                     Rt(ix,iy,iz)=Rt(ix,iy,iz)+s*R;
%                     
%                     % end;
% 
%                 end
%             end
%         end
%     end
%     
% end;

%figure(2);scatter3(Det_hit(1,:),Det_hit(2,:),Det_hit(3,:));title('scatter of detector hit');xlabel('x');ylabel('y');zlabel('z');

S2=flipdim(S,2);
S3=(S+S2)*9/10;
Rt_weight=(Rt+flipdim(Rt,2))./(S+S2);

%Process z sequence
Stot=zeros(size(S3));
for i=1:n_seq_z
    [act_k,mat_k]=get_indices(nz,nz,round(seq_z(i)/dz));
    Stot(:,:,act_k)=Stot(:,:,act_k)+S3(:,:,mat_k);
end;
Stot=Stot/n_seq_x/n_seq_z;
%view_nii(make_nii(Stot*100))
sum(Stot(:));


width=150;
chesttonipple=110;
thickness=55;


%create boolean mask and simultanously coordinate tables
breast=false(nx,ny,nz);
for i=1:nx
    for j=1:ny
        for k=1:nz
            if ((x(i)/(width/2))^2+((685-z(k))/chesttonipple)^2)<=1 && z(k)<675 && y(j)>-thickness/2&& y(j)<thickness/2
                breast(i,j,k)=true;
            end
        end
    end
end

msk=(abs(Y)<27.5)&breast;%create mask
Sens=[mean(Stot(msk)) mean(Stot((abs(Y)<15)&(abs(Z)<650)&(abs(Z)>620)&(abs(X)<15)))  mean(Rt_weight(msk&~isnan(Rt_weight))) mean(Rt_weight((abs(Y)<15)&(abs(Z)<650)&(abs(Z)>620)&(abs(X)<15)&~isnan(Rt_weight)))];
%toc;