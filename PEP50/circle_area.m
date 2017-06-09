%TO DO:
%-scale the width of the pinhole depending on the angles
%-find a better optimization
%-check the grid
%- HOW THE FUCK DOES THE GRID EVEN WORK
%- maybe vectorize



function sphere_area = circle_area(hit_pin_xyz, grid_col, grid_row)
    %we nemen de unitsphere als een bol met straal 1mm.
    %en xyz ook in mm. de eerste colom van hit_pin bevat allen de xyz van de geraakte voxel
    
    %determine the xyz coordinates of the voxel being hit
    v = hit_pin_xyz(1:3,1);
    hit_pin = hit_pin_xyz(:,2:end);
    sz = size(hit_pin,2);
    for i=1:sz
        hit_pin(1:3,i) = hit_pin(1:3,i) - v(1:3); 
    end
    
    %here we calculate phi and theta note, phi is the angle with (x,y)
    phi = atan( hit_pin(2,:) ./ hit_pin(1,:) ) + (hit_pin(1,:) < 0) .* (hit_pin(2,:) > 0)*pi/2 - (hit_pin(1,:) < 0) .* (hit_pin(2,:) < 0)*pi/2;
    theta = atan( hit_pin(3,:) ./ sqrt( hit_pin(1,:).^2 + hit_pin(2,:).^2) );
    
    %now we want to calculate the area on the circle, we do this by
    %calculating the borders. we take 12 as a estimate to make sure we hit
    %all the grids.
    
    %determine the width of the pinholes
    r_pin = hit_pin(6,:);
    %but now in angles.
    r_voxel = 0.1;
    r_pin = atan(r_voxel + r_pin ./ sqrt( hit_pin(1,:).^2 + hit_pin(2,:).^2 + hit_pin(3,:).^2));
    
    %maak het aantal stappen moet overlegd worden
    steps = max(grid_row,grid_col)*4;
    gamma = linspace(2*pi/steps, 2*pi, steps);
    
    
    %determine the angles of the surface vector
    pin_phi = hit_pin(4,:);
    pin_theta = hit_pin(5,:);
    
    %begin making the area that is been seen by  the pinhole
    area_border = zeros(sz, steps, 2);
    for i=1:sz
        for j=1:steps
            area_border(i,j,:) = [(phi(i) + cos(gamma(j)).*r_pin(i)) (theta(i) + sin(gamma(j)).*r_pin(i))];
        end
    end
    
    %begin making a grid area of the area border.
    
    %row = niet x2
    %grid phi = row
    area_border(:,:,1) = round(area_border(:,:,1) / (2 * pi) * grid_row)/grid_row*2*pi;
    area_border(:,:,2) = round(area_border(:,:,2) / pi * grid_col)/grid_col*pi;
    sphere_area = area_border;
    %%{
    figure;
    hold on;
    for i=1:sz
        plot(area_border(i,:,1),area_border(i,:,2));
    end
    axis([-pi pi -pi/2 pi/2]);
    grid on;
    %%}
end
