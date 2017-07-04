function [transpin] = all_pinhole_positions(M,R_transaxial,N_Pos_rev,N_Pos_tot,scanning_input)
%M = [Lx Ly Lz; volume of spiral 
%R_transaxial; raidus of spiral
%N_Pos_rev; number of steps per revolution of spiral
%N_Pos_tot; total number of steps
%scanning_input is the entire pinhole information

phantom_size = [M(1) M(2) M(3)];
Z_step = (phantom_size(3))/(N_Pos_tot-1);

X=[];Y=[];Z=-phantom_size(3)/2;
for l=1:N_Pos_tot
    X = [X,R_transaxial*cos(2*pi/N_Pos_rev*(l-1))];
    Y = [Y,R_transaxial*sin(2*pi/N_Pos_rev*(l-1))];
    Z = [Z,Z(end)+Z_step];
end
Z = Z(1:end-1);

for i=1:size(Z,2)
    if abs(Z(i))>abs(phantom_size(3)/2)-3 && Z(i)<0
        Z(i)=-phantom_size(3)/2+3;
    elseif abs(Z(i))>abs(phantom_size(3)/2)-3 && Z(i)>0
        Z(i)=phantom_size(3)/2-3;
    end
end
scannerPositions = vertcat(X, Y, Z)';

transpin = zeros(75*length(scannerPositions),7);

pinholes = scanning_input;
for i = 1:length(scannerPositions);
    
    transpin(75*(i-1)+1:75*i,1) = pinholes(:,1) + scannerPositions(i,1);
    transpin(75*(i-1)+1:75*i,2) = pinholes(:,2) + scannerPositions(i,2);
    transpin(75*(i-1)+1:75*i,3) = pinholes(:,3) + scannerPositions(i,3);
    transpin(75*(i-1)+1:75*i,4:7) = pinholes(:,4:7) ;  % draws the pinholes

end
end 
