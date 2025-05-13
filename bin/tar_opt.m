function [ R,mod ] = tar_opt( p,xf,datf,anz,plen,we,lenx,lenxf,actx,actx_w,ntol)

%disp(['plenanz: ' num2str(plen) ' ' num2str(anz) ' ' num2str(length(p))])
P=reshape(p,plen,anz);


modf=zeros(size(xf));
for i=1:anz
    modf=modf+mod_ga(xf,P(:,i));
end


% M=zeros(length(xf),anz);
% for i=1:anz
%     M(:,i)=mod_ga(xf,P(:,i));
% end
% modf=M*ones(anz,1);

mod     = modf(actx);
dat     = datf(actx);
% if we(2)>0
    modw    = modf(actx_w);
    datw    = datf(actx_w);
% end

% coeff=1/(max(abs(dat)));
%Modellanpassung
mdiff=dat-mod;
if ntol>0
    id=abs(mdiff)<=ntol;            
            
    mdiff(id)=0;
end

% R1=we(1)*coeff*(mdiff);
R1=we(1)*(mdiff(:));

%Modell kleiner als Daten in Fenster
% R2=we(2)*0.1*coeff*min(datw-modw,0);
% if we(2)>0
R2=we(2)*0.1*min(datw-modw,0);
R=[R1(:); R2(:)];
% end

%Modell kleiner als Daten in grossem Fenster
% if lenxf-lenx>0
%     R3=we(2)*0.1*coeff*(lenx/(lenxf-lenx))*min(datcomp-modcomp,0);
% end
%Symmetrie approximativ Modell 
%R3=50*max(abs(P(2,:)-P(3,:))-0.5,0);

% if lenxf-lenx>0
%     R=[R1(:); R2(:); R3(:)];%; R3(:)];
% else
%     R=[R1(:); R2(:)];
% end
end

