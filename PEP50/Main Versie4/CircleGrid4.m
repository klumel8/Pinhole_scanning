function [gridcir]=CircleGrid4(gridrow,gridcol,part)


rowleng=1.*pi/gridrow;
colleng=2.*pi/gridcol;


numberSteps = max(gridrow,gridcol)*15;
dg = 2*pi/numberSteps;
numberCircles = gridcol*((gridrow-1)/2+1);
gridcir=false(gridcol,gridrow,floor(numberCircles/part(2))/100); %devide by 100 for large res
a = floor(numberCircles*(part(1)-1)/part(2))+1;
start = a-1;
final = floor(numberCircles*part(1)/part(2))

%h = waitbar(0,'Creating circles...');

i = 1:gridcol;
pm = pi/gridcol + (i-1)*2*pi/gridcol;

j = 0:(gridrow-1)/2;
tm = 1/2*pi/gridrow + j*pi/gridrow;
        
All = meshgrid(pm,tm);
        

for k = a:floor((final-a))/100+a %devide by 100
    %weight(a)=cos(tm);%(2*pi./gridcol*(cos(tm-pi/(gridrow-1)/2)-cos(tm+pi/(gridrow-1)/2))/(4*pi));
    %{
    if mod(k,1000) == 0
        disp(k)
    end
    %}
    for gamma = 1:numberSteps
        x = (cos(dg*gamma)*sin(All(k))*cos(All(k))-sin(dg*gamma)*sin(All(k)));
        y=(cos(dg*gamma)*sin(All(k))*sin(All(k))+sin(dg*gamma)*cos(All(k)));
        z=(cos(dg*gamma)*cos(All(k)));
        r=sqrt((x).^2+(y).^2+(z).^2);
        thetal=acos(z./r);
        phil=atan2(y,x);
        
        if phil<0
            phil=phil+2*pi;
        end
        
        gridcir(ceil(phil/colleng),ceil(thetal/rowleng),k-start)=1;
        %gridcir(ceil(thetal/rowleng),ceil(phil/colleng),i)+1;
        
    end
    % gridcir(:,:,i)=gridcir(:,:,i)/sum(sum(gridcir(:,:,i)));
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
