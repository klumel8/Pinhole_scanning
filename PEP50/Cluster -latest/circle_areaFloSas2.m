%TO DO:

%-scale the width of the pinhole depending on the angles
%-find a better optimization
%-check the grid
%- HOW THE FUCK DOES THE GRID EVEN WORK
%- maybe vectorize


function sphere_area = circle_area(hit_pin_xyz, grid_row, grid_col)
    %close all
    %we nemen de unitsphere als een bol met straal 1mm.
    %en xyz ook in mm. de eerste colom van hit_pin bevat allen de xyz van de geraakte voxel

    %determine the xyz coordinates of the voxel being hit (first colomn of
    %the pinhole list)
    v = hit_pin_xyz(1:3,1);
    %remove the voxel coordinates column from the matrix.
    hit_pin = hit_pin_xyz(:,2:end);
    sz = size(hit_pin,2);
    for i=1:sz
        hit_pin(1:3,i) = -(hit_pin(1:3,i) - v(1:3)); 
    end
    
    %here we calculate phi and theta note, phi is the angle with (x,y)
    phiC = atan2(hit_pin(2,:), hit_pin(1,:))+pi;
    %phi = phi .* (phi >= 0) + (phi + 2 * pi) .* (phi < 0); 
    thetaC = abs(atan( hit_pin(3,:) ./ sqrt( hit_pin(1,:).^2 + hit_pin(2,:).^2))-pi/2);
    
    %now we want to calculate the area on the circle, we do this by
    %calculating the borders. we take 12 as a estimate to make sure we hit
    %all the grids.
    %determine the angles of the surface vector
    pin_phi = hit_pin(4,:)';
    pin_theta = hit_pin(5,:)';
    r_pin = hit_pin(6,:)';
    lengte = sqrt([hit_pin(1,:)'.^2 + hit_pin(2,:)'.^2 + hit_pin(3,:)'.^2]);
    mess1 = -1.*[hit_pin(1,:)'./lengte , hit_pin(2,:)'./lengte, hit_pin(3,:)'./lengte];
    mess2 = [cos(pin_phi).*sin(pin_theta), sin(pin_phi).*sin(pin_theta), cos(pin_theta)]; 
    angle = acos((mess1(:,1).*mess2(:,1)+mess1(:,2).*mess2(:,2)+mess1(:,3).*mess2(:,3)));
    r_pin =(r_pin.*abs(sin(pi/2 - angle))./(lengte + r_pin.*abs(cos(pi/2-angle))));
    r_pin = r_pin';
    r_pin = atan(r_pin);
    
    %maak het aantal stappen moet overlegd worden
    steps = 500;
    h_step = round(steps/2);
    gamma = linspace(2*pi/steps, 2*pi, h_step);
    
    %begin making the area that is been seen by  the pinhole
    area_border = zeros(h_step, sz, 2); %steps
    [Phi,Gamma] = meshgrid(phiC,gamma);
    [Theta,~] = meshgrid(thetaC,gamma);
    [R_pin,~] = meshgrid(r_pin,gamma);
    
    phiC = phiC-pi;
    thetaC = -(thetaC+pi/2);
    
    Phim = Phi-pi;
    Thetam = -(Theta + pi/2);
    pos = ones(h_step,sz,3);                          
        
    pos(1:h_step,:,1) = R_pin.*((cos(Gamma).*sin(Thetam).*cos(Phim)-sin(Gamma).*sin(Phim)))+cos(Phi).*sin(Theta);
    pos(1:h_step,:,2) = R_pin.*(cos(Gamma).*sin(Thetam).*sin(Phim)+sin(Gamma).*cos(Phim))+sin(Phi).*sin(Theta);
    pos(1:h_step,:,3) = R_pin.*cos(Gamma).*cos(Thetam)+cos(Theta);
    
    %pos(1:h_step,:,1) = pos(1:h_step,:,1)+cos(Phi).*sin(Theta);                      
    %pos(1:h_step,:,2) = pos(1:h_step,:,2)+sin(Phi).*sin(Theta);
    %pos(1:h_step,:,3) = pos(1:h_step,:,3)+cos(Theta);
    
    %plot3(pos(1:h_step,:,1),pos(1:h_step,:,2),pos(1:h_step,3));
    
    %area_border(:,:,2) = atan2(pos(:,:,2), pos(:,:,1))+ pi;
    %phi = phi .* (phi >= 0) + (phi + 2 * pi) .* (phi < 0); 
    %area_border(:,:,1) = abs(atan( pos(:,:,3) ./ sqrt( pos(:,:,1).^2 + pos(:,:,2).^2))-pi/2);
    
    area_border(:,:,2) = floor((atan2(pos(:,:,2), pos(:,:,1))+ pi)/2/pi*grid_col + 1); %Phi?
    area_border(:,:,1) = floor(abs(atan( pos(:,:,3) ./ sqrt( pos(:,:,1).^2 + pos(:,:,2).^2))-pi/2)/pi*grid_row + 1); %Theta?
    %area_border(:,:,2) = (atan2(pos(:,:,2), pos(:,:,1))+ pi); %Phi?
    %area_border(:,:,1) = abs(atan( pos(:,:,3) ./ sqrt( pos(:,:,1).^2 + pos(:,:,2).^2))-pi/2); %Theta?
    
    
    %area_border(:,:,2) = floor(area_border(:,:,2)/2/pi*grid_col + 1); %Phi?
    %area_border(:,:,1) = floor(area_border(:,:,1)/pi*grid_row + 1); %Theta?
    
    %plot(area_border(1:h_step,:,2),area_border(1:h_step,:,1),'.');
    
    %{
    for i=1:sz
        for j=1:steps
            area_border(i,j,:) = [(phi(i) + cos(gamma(j)).*r_pin(i)) (theta(i) + sin(gamma(j)).*r_pin(i))];
        end
    end
    area_border - area_border1
    %}
    
    %begin making a grid area of the area border.
    
    %row = niet x2
    %grid phi = row
    %figure
    %plot(area_border(1:h_step,:,2),area_border(1:h_step,:,1),'.');
    
    
    %{
    hold on;
    for i=1:sz
        plot(area_border(:,i,2),area_border(:,i,1), '.');
        xlabel('phi')
        ylabel('theta')
    end
    axis([0 2*pi 0 pi]);
    grid on;
    %}
    
    %now make the border area into a grid matrix;
    border_grid = false(grid_row,grid_col);
    %size(border_grid)
    for i=1:sz
        %for j=1:h_step %steps
            %r_phi = floor(area_border(j,i,2)/2/pi*grid_col + 1);
            %r_theta = floor(area_border(j,i,1)/pi*grid_row + 1);
            border_grid(area_border(:,i,1),area_border(:,i,2)) = true;
        %end
    end
    %size(border_grid)
    
%    colormap('gray');
plot2 = border_grid(:,:);
%imagesc(plot2);

    
    sphere_area = border_grid(1:grid_row,1:grid_col);
end