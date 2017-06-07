%Hey
%creating cones section
%Matrix M (position vector cone, angle vector cone, maximum angle)
M = [1 1 1 -1 -1 -1 0.05*pi; -1 -1 1 1 1 -1 0.05*pi];
Z = ['none'];
for i = 1:size(M,1)
    len = sqrt((M(i,4)-M(i,1))^2 +(M(i,5)-M(i,2))^2 + (M(i,6) - M(i,3))^2);
    Cone([M(i,1) M(i,2) M(i,3)],[M(i,4) M(i,5) M(i,6)],[0 len*sin(M(i,7))],60,Z,0,1);
    hold all
end 
%% 
%creating points
N = 1000;
x = -1 + 2*rand(N,1);
y = -1 + 2*rand(N,1);
z = -1 + 2*rand(N,1);

%% 
%checking if point is inside cone(s)
int = ones(N,1);
in =ones(N,size(M,1));
for j = 1:size(M,1)
x0 = x - M(j,1);
y0 = y - M(j,2);
z0 = z - M(j,3);
v = [M(j,4) M(j,5) M(j,6)];
k = norm(v);
for i = 1:length(x0);
in(i,j) = (M(j,7) > acos( dot([x0(i) y0(i) z0(i)],v)/(sqrt(x0(i)^2 + y0(i)^2 + z0(i)^2)*k)));
end 
int = int.*in(:,j);
outt = ~logical(int);
end 
%% 
%plotting points
plot3(x(logical(int)),y(logical(int)),z(logical(int)),'or', x(logical(outt)),y(logical(outt)),z(logical(outt)),'og');
xlabel('x'); ylabel('y'); zlabel('z');
grid on;