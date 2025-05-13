function [res] = job_fcn(j,p,tg,noiselvl,ClusterIdx,fev_tol,info,init_idx,calcWin)
res.success=0; %0=failed, 1=success
res.ClusterIdx=ClusterIdx;
res.init_idx=init_idx;
res.calcWin=calcWin;



%% tar fun details
% % tar_opt_par(p,...                          %para
% %     j.timegrid_shi_red,...                 %shifted time grid reduced to wide window
% %     j.res_red_scaled,...                   %signal scaled to max 1 reduced to wide window
% %     [p.idx_difStartInNeg p.idx_difLen],... % narrow window information for determination in wide window
% %     [1 6],...                              %parameter size info
% %     we,...                                 %weights
% %     ntol);                                 %comp-wise absolute error tolerance

%% optimizations
p0=j.p0;

if j.lb(2)>=1E-5 && j.lb(3)>=1E-5
    if 1
        we=[1 0];
        %Zielfunktion in gewuenschte Form bringen (Abhaengigkeit von nur einer Variable)
        fun2=@(pa) tar_opt_par(pa,...
            j.timegrid_shi_red,...
            j.res_red_scaled,...
            [p.idx_difStartInNeg p.idx_difLen],...
            [6 1],...
            we,...
            0);
        options=optimoptions(@lsqnonlin,...
            'Display','off',...
            'MaxIter',100,...
            'MaxFunEvals',10000,...
            'TolFun',1e-10,...
            'TolX',1e-10);
        pOpt=lsqnonlin(fun2,p0,j.lb,j.ub,options);
    else
        pOpt=p0;
    end
    
    we=[1 10];
    fun2=@(pa) tar_opt_par(pa,...
        j.timegrid_shi_red,...
        j.res_red_scaled,...
        [p.idx_difStartInNeg p.idx_difLen],...
        [6 1],...
        we,...
        0);
    options=optimoptions(@lsqnonlin,...
        'Display','off',...
        'MaxIter',100,...
        'MaxFunEvals',10000,...
        'FinDiffType','central',...
        'TolFun',1e-10,...
        'TolX',1e-10);
    pOpt=lsqnonlin(fun2,pOpt,j.lb,j.ub,options);
else
    pOpt=p0;
end

%% eval target fcn
we=[1 1];
fvec=tar_opt_par(pOpt,...
    j.timegrid_shi_red,...
    j.res_red_scaled,...
    [p.idx_difStartInNeg p.idx_difLen],...
    [6 1],...
    we,...
    noiselvl/j.res_red_scal);

fev=sum(fvec.^2);


%% debug
debug=0;
if debug && ClusterIdx==1
    figure(555);
    clf
%     subplot(2,1,1)
    
    plot(j.timegrid_shi_red,j.res_red_scaled,'b');
    hold on
    plot(j.timegrid_shi_red([1 end]), noiselvl/j.res_red_scal*[1 1],'r--')
    mod=mod_ga(j.timegrid_shi_red,pOpt);
    mod0=mod_ga(j.timegrid_shi_red,p0);
    plot(j.timegrid_shi_red,mod0,'k')
    plot(j.timegrid_shi_red,mod,'r')
    title(['par ' num2str(fev) ' of ' num2str(fev_tol) ' at ' num2str(p.center)])
    1;
end



%% postprocess
%de-shift
pOpt(1)=pOpt(1)+p.center;
p0(1)=p0(1)+p.center;
%de-scale
pOpt(4)=j.res_red_scal*pOpt(4);
p0(4)=j.res_red_scal*p0(4);


res.mod=mod_ga(tg,pOpt); %unscaled model on tg
res.pOpt=pOpt;
res.p0=p0;
res.fev=fev; %relative residuals

%% check if fit is ok

if pOpt(4)>noiselvl && fev<fev_tol
    if info>2
        disp(['<> ADD peak at ' num2str(pOpt(1),4)]);
    end
    
    res.success=1;
else
    if info>2
        disp(['<> IGNORE fitted peak: hi= ' num2str(pOpt(4)) '<=' num2str((noiselvl)) '; fev= ' num2str(fev) '>=' num2str(fev_tol)])
    end
    
    %fcn returns with res.success=0;
end






end




















