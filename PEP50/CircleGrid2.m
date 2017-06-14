%{
function gridcir=CircleGrid(Cp,gridrow,gridcol)


rowleng=1.*pi/gridrow;
colleng=2.*pi/gridcol;

c=size(Cp);
gridcir=zeros(gridrow,gridcol,c(3));

for i=1:c(3) 
    for j=1:c(1)
             
        thetal=Cp(j,1,i);
        
        phil=Cp(j,2,i);
        if phil<0
            phil=phil+2*pi;
        end;                  
               
        gridcir(floor(thetal/rowleng)+1,floor(phil/colleng)+1,i)=1;
        
    end
  
    %a=(gridrow+1)/2;
   % b=c(3)+1;
    
    %gridcir(a,:,b)=1;
    
end
%}

function [gridcir]=CircleGrid(gridrow,gridcol)


rowleng=1.*pi/gridrow;
colleng=2.*pi/gridcol;

a = 0;
Nc = max(gridrow,gridcol)*5;
dg = 2*pi/Nc;
ngrid = gridcol*((gridrow-1)/2+1);
gridcir=false(ngrid,gridrow,gridcol);



%h = waitbar(0,'Creating circles...');
weight=zeros(ngrid,1);
for i = 1:gridcol
    pm = pi/gridcol + (i-1)*2*pi/gridcol;
 %   waitbar((i-1)/(gridcol),h);
    % pm/pi*180
    for j = 0:(gridrow-1)/2-1;
        tm = 1/2*pi/gridrow + j*pi/gridrow;
        
        a = a+1;
        weight(a)=cos(tm);%(2*pi./gridcol*(cos(tm-pi/(gridrow-1)/2)-cos(tm+pi/(gridrow-1)/2))/(4*pi));
        
        for gamma = 1:Nc
            x = (cos(dg*gamma)*sin(tm)*cos(pm)-sin(dg*gamma)*sin(pm));
            y=(cos(dg*gamma)*sin(tm)*sin(pm)+sin(dg*gamma)*cos(pm));
            z=(cos(dg*gamma)*cos(tm));
            r=sqrt((x).^2+(y).^2+(z).^2);
            thetal=acos(z./r);
            phil=atan2(y,x);
            
            if phil<0
                phil=phil+2*pi;
            end;
            
            gridcir(a,ceil(thetal/rowleng),ceil(phil/colleng))=1;%gridcir(ceil(thetal/rowleng),ceil(phil/colleng),i)+1;
        end
    end    
    % gridcir(:,:,i)=gridcir(:,:,i)/sum(sum(gridcir(:,:,i)));
end
 gridcir = permute(gridcir,[3,2,1]);
end
%close(h);
%weight=weight/sum(weight);
%no equatior, is not necessary. already covered by the other squares so
%long as grid is fine enough.
%
%a=(gridrow+1)/2;
%b=c(3)+1;
%gridcir(a,:,b)=1;
%gridcir(:,:,b)=gridcir(:,:,b)/sum(sum(gridcir(:,:,b)));
