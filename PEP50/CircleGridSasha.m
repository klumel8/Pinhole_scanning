%TO DO:
%-scale the width of the pinhole depending on the angles
%-find a better optimization
%-check the grid
%- HOW THE FUCK DOES THE GRID EVEN WORK
%- maybe vectorize



function sphere_area = circle_area(hit_pin_xyz, grid_col, grid_row)
    %we nemen de unitsphere als een bol met straal 1mm.
    %en xyz ook in mm. de eerste colom van hit_pin bevat allen de xyz van de geraakte voxel
    
    %determine the xyz coordinates of the voxel being hit (first colomn of
    %the pinhole list)
    v = hit_pin_xyz(1:3,1);
    %remove the voxel coordinates column from the matrix.
    hit_pin = hit_pin_xyz(:,2:end);
    sz = size(hit_pin,2);
    for i=1:sz
        hit_pin(1:3,i) = hit_pin(1:3,i) - v(1:3); 
    end
    
    %here we calculate phi and theta note, phi is the angle with (x,y)
    phi = atan( hit_pin(2,:) ./ hit_pin(1,:) ) + (hit_pin(1,:) < 0) * pi + pi/2;
    theta = abs(atan( hit_pin(3,:) ./ sqrt( hit_pin(1,:).^2 + hit_pin(2,:).^2))-pi/2);
    
    %now we want to calculate the area on the circle, we do this by
    %calculating the borders. we take 12 as a estimate to make sure we hit
    %all the grids.
    %determine the angles of the surface vector
    pin_phi = hit_pin(4,:);
    pin_theta = hit_pin(5,:);
    
    %do somethign artur told me
    lengte = norm([hit_pin(1,:)  hit_pin(2,:) hit_pin(3,:)]);
    mess1 = -1 * [hit_pin(1,:)  hit_pin(2,:) hit_pin(3,:)] / lengte;
    mess2 = [cos(pin_phi).*sin(pin_theta) sin(pin_phi).*sin(pin_theta) cos(pin_theta)];
    delta = pi/2 - acos(mess1(:,1) .* mess2(:,1) + mess1(:,2) .* mess2(:,2) + mess1(:,3) .* mess2(:,3));
    
    %determine the width of the pinholes
    r_pin = hit_pin(6,:)*5;
    r_pin = (r_pin * lengte) ./ (lengte + r_pin.*cos(delta));
    %but now in angles.
    r_voxel = 0.0;
    r_pin = atan(r_voxel + r_pin ./ lengte);
    
    %determine the width of the pinholes
    r_pin = hit_pin(6,:);
    %but now in angles.
    r_voxel = 0.0;
    r_pin = atan(r_voxel + r_pin ./ sqrt( hit_pin(1,:).^2 + hit_pin(2,:).^2 + hit_pin(3,:).^2));
    
    %maak het aantal stappen moet overlegd worden
    steps = 100;
    h_step = round(steps/2);
    gamma = linspace(2*pi/steps, 2*pi, steps);
    
    %begin making the area that is been seen by  the pinhole
    area_border = zeros(steps, sz, 2);
    [X,Y] = meshgrid(1:sz,1:h_step);
    area_border(1:h_step,:,2) = phi(X) + cos(gamma(Y)).*r_pin(X);
    area_border(1:h_step,:,1) = theta(X) + sin(gamma(Y)).*r_pin(X);
    area_border(h_step+1:end,:,2) = phi(X) + 0.8*cos(gamma(Y)).*r_pin(X);
    area_border(h_step+1:end,:,1) = theta(X) + 0.8*sin(gamma(Y)).*r_pin(X);
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
    area_border(:,:,1) = round(area_border(:,:,1) / pi * grid_row)/grid_row*pi;
    area_border(:,:,2) = round(area_border(:,:,2) / 2 / pi * grid_col)/grid_col*2*pi;
    %%{
    hold on;
    for i=1:sz
        plot(area_border(:,i,2),area_border(:,i,1));
    end
    axis([0 2*pi 0 pi]);
    grid on;
    %}
    %now make the border area into a grid matrix;
    border_grid = false(grid_row,grid_col);
    %border_grid1 = border_grid;
    %r_phi = ceil(area_border(1:steps,1:sz,2)/2/pi*grid_col);
    %r_theta = ceil(area_border(1:steps,1:sz,1)/pi*grid_row);
    %r_theta = reshape(r_theta,1,size(r_theta,1)*size(r_theta,2));
    %r_phi = reshape(r_phi,1,size(r_phi,1)*size(r_phi,2));
    %ray = [r_theta; r_phi];
    %border_grid1(ray) = true;
    for i=1:sz
        for j=1:steps
            r_phi = ceil(area_border(j,i,2)/2/pi*grid_col);
            r_theta = ceil(area_border(j,i,1)/pi*grid_row);
            border_grid(r_theta,r_phi) = true;
        end
    end
    %isequal(border_grid,border_grid1)
    sphere_area = border_grid;
end

%d = (pi/2) -acos(-[
