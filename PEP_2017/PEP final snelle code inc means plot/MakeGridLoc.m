function gridloc=MakeGridLoc(gridrow,gridcol,phivoxvox,thetavox,FOV,Npin)

fov=length(FOV);
phivoxg=(180./pi).*phivoxvox;
thetag=(180./pi).*thetavox;
gridloc=zeros(fov,3);
gridloc(:,3)=FOV;

rowleng=180./gridrow;
colleng=360./gridcol;


gridloc(:,1)=floor(thetag/rowleng)+1;
gridloc(:,2)=floor(phivoxg/colleng)+1;

