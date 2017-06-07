function mask = makemask(nvoxX,nvoxY,nvoxZ,voxelsize,radius_mouse)

mask = zeros(nvoxX,nvoxY,nvoxZ);
center = size(mask)/2+0.5;
radius_mouse = radius_mouse/voxelsize;

for x=1:size(mask,1)
for y=1:size(mask,2)
for z=1:size(mask,3)

    if sqrt((x-center(1))^2+(y-center(2))^2)<=radius_mouse
        mask(x,y,z)=1; 
    end
            
end
end
end

end