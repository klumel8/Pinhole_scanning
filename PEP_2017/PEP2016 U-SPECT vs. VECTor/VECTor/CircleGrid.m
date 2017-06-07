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
  
    a=(gridrow+1)/2;
    b=c(3)+1;
    
    gridcir(a,:,b)=1;
    
end

