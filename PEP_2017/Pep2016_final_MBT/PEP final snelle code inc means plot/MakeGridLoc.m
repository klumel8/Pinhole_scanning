function gridloc=MakeGridLoc(gridrow,gridcol,phivox,thetavox)
phivoxg=(180./pi).*phivox;
thetag=(180./pi).*thetavox;
gridloc=zeros(length(phivox),2); 

rowleng=180./gridrow;
colleng=360./gridcol;

gridloc(:,1)=floor(thetag/rowleng)+1;
gridloc(:,2)=floor(phivoxg/colleng)+1;
