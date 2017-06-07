 pp1=gridcir2(1,:,:);
 pp2=squeeze(pp1); % 19x32
 size(pp2)
 
 x=1:32;
 y=1:19;
 
 xx=zeros(19,32);
 
 for i=1:32
     xx(1:length(x(pp2(i,:))),i)=x(pp2(:,i));
 end
     
%      for j=1:19
     