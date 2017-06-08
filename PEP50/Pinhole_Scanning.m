%Hey
%creating cones section
%Matrix M (position vector cone, angle vector cone, maximum angle)
%for x y z phi theta radius alpha
M = [1 1 1 1.25*pi 0.75*pi 0.5 0.05*pi];
C= zeros(size(M));
C(:,1:3) = M(:,1:3) - (M(:,6)./tan(M(:,7))).*[cos(M(:,4)).*sin(M(:,5)) sin(M(:,4)).*sin(M(:,5)) cos(M(:,5))];
C(:,4:end) = M(:,4:end);
Z = ['none'];
for i = 1:size(M,1)
    Cone([C(i,1) C(i,2) C(i,3)],[M(i,1) M(i,2) M(i,3)],[0 0.9],60,Z,0,1);
    hold all
end 
%% 
%creating points
N = 1000;
x = -5 + 5*rand(N,1);
y = -5 + 5*rand(N,1);
z = -5 + 5*rand(N,1);

%% 
%checking if point is inside cone(s)
int = ones(N,1);
in =ones(N,size(C,1));
for j = 1:size(C,1)
x0 = x - C(j,1);
y0 = y - C(j,2);
z0 = z - C(j,3);
v = [cos(C(:,4)).*sin(C(:,5)) sin(C(:,4)).*sin(C(:,5)) cos(C(:,5))]';
in(:,j) = (C(j,7) > acos( [x0 y0 z0]*v./(sqrt(sum(abs([x0 y0 z0]).^2,2))))); 
int = int.*in(:,j);
outt = ~logical(int);
end 
%% 
%plotting points
plot3(x(logical(int)),y(logical(int)),z(logical(int)),'or', x(logical(outt)),y(logical(outt)),z(logical(outt)),'og');
xlabel('x'); ylabel('y'); zlabel('z');
hold all
quiver3(C(:,1),C(:,2),C(:,3),10*cos(M(:,4)).*sin(M(:,5)), 10*sin(M(:,4)).*sin(M(:,5)), 10*cos(M(:,5)))
grid on;