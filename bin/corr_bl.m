function [RESBL] = corr_bl(gx,gc,ng,d,N)
scal=max(abs(gc));
gcs=gc./scal;
RESBL.scal=scal;

%%
mod=zeros(size(gx));

rerun=2;
M=zeros(rerun,length(gx));
if ng>0
    centeridx  = round(linspace(1,length(gc),ng));
    dx=abs(gx(centeridx(2))-gx(centeridx(1)));


    for i=1:rerun
        center    = linspace(gx(1),gx(end),ng)+(i-1)/(rerun)*dx;
        hi        = 0.1*ones(ng,1);
        if center(end)>gx(end)
            center(end)=gx(end);
        end
        if center(1)>gx(1)
            center = [gx(1) center];
            hi     = [hi(1); hi];
        end
        
        fun=@(para_hi) tar_bl( para_hi,gx,sgfilt(d,N,gcs),center);
        options = optimoptions(@lsqnonlin,...
            'TolFun',1E-8,...
            'TolX',1E-8,...
            'MaxFunEvals',2000,...
            'MaxIter',50,...
            'Display','iter');
        hi_opt=lsqnonlin(fun,hi,zeros(size(hi)),[],options);
        curmod=spline(center,hi_opt,gx);
        M(i,:)=curmod;
    end

    gcs=gcs-mean(M)';

RESBL.bl.mod=mod*scal;
RESBL.gcs=gc-mod*scal;

end

