function C = CircleGen(gridrow,gridcol)
a = 0;
Nc = 4*max(gridrow,gridcol);
dg = 2*pi/Nc;
ngrid = gridcol*(gridrow-1)/2;
C = zeros((gridrow+1)/2,3,ngrid);

for i = 1:gridcol
    pm = pi/gridcol + i*2*pi/gridcol;
    for j = 0:(gridrow-1)/2-1;
        tm = 1/2*pi/gridrow + j*pi/gridrow;
        a = a+1;
        for gamma = 1:Nc
            C(gamma,:,a) = [(cos(dg*gamma)*sin(tm)*cos(pm)-sin(dg*gamma)*sin(pm)) (cos(dg*gamma)*sin(tm)*sin(pm)+sin(dg*gamma)*cos(pm)) (cos(dg*gamma)*cos(tm))];
        end
    end
end
    
