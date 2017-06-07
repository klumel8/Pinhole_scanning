function Cp=Polariseer(C)
c=size(C);

for i=1:c(3)
    for j=1:c(1)
        x=C(j,1,i);
        y=C(j,2,i);
        z=C(j,3,i);
        r=sqrt((x).^2+(y).^2+(z).^2);
        
        Cp(j,1,i)=acos(z./r);
        Cp(j,2,i)=atan2(y,x);
    end
end