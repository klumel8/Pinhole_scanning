% create the new changed pinhole coordinates of the whole translated
% system
k = 4;  %number_of_positions

kut = linspace(0,pi, k);
pinholes_rot = pinholes(:,1:3)'; %x, y,z coordinaten voor het roteren

for j = 2:k;  
    i = kut(k); 
    rot = [cos(i) sin(i) 0; -sin(i) cos(i) 0 ; 0 0 1]; % rotatie matrix
  
    A = (rot*pinholes_rot)'; 
    B = pinholes(1:75, 4:7);
    B(:,4) = B(:, 4 )+ i;  % hoek verplaatsing
    A = horzcat( A, B);
    pinholes = vertcat(pinholes, A);
end
clear A B i j k kut rot