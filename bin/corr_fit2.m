function [RES] = corr_fit2(para,gx,gcs,d,N,noiseLVL)


%% set dynamic parameters
% info level for outputs: 0=none, 1=partially, 2=all
info=para.info;
% maximal peak number
maxpeak_nr=para.maxpeak_nr;
% abort error value
errTol=para.errTol;
% savitzky golay parameter
epsi=para.epsi;
delta=para.delta;
% lack of fit tolerance to accept a peak 
fev_tol=para.fev_tol;
% number of runs to check for overlapping peaks
overlap_runs=para.overlap_runs;
% amount of overlap for detecting peak intersections (0..1)
p_correl=para.peakoverlap;

%% Set constant parameters
p_area=3;


%% START
RES.peaks=[];

% scaling
scal=max(abs(gcs));
gcs=gcs./scal;
RES.scal=scal;

% model init
gcmod=zeros(length(gx),1);

% diff init
err=(trapz(gcs)-trapz(gcmod))/trapz(gcs);




%% FIT
num=0;
if info>1
    disp(['<> ' num2str(num) '; error: ' num2str(err) '; errortol: ' num2str(errTol)])
end
nocurve_prev=0;
% active blocks (will be deactivated)
idxBlock=true(size(gcs));
idxfull=1:length(gcs);


while err>errTol && scal*max(gcs-gcmod)>3*noiseLVL
    
    
    num=num+1;
    % max peak number reached?
    if length(RES.peaks)>maxpeak_nr-1
        break
    end
    
    % current residual profile
    currGC=gcs-gcmod;
    % Savitzky Golay smoothing
    [d_smooth,~,~,~]= sgfilt(d,N,currGC);
    scal1=max((currGC));
    % find max difference
    [Vymin,idxtemp]=max(d_smooth(idxBlock)./scal1);
    idxreduce=idxfull(idxBlock);
    Vxminidx=idxreduce(idxtemp);
    Vxmin=gx(Vxminidx);
    
    
    ign=0;
    curve_added=0;
    % abort if para.retry subsequent trys no peak were found
    if nocurve_prev==para.retry
        break
    end
    
    if ~isempty(Vxmin)
        currGC=(gcs-gcmod);
        scal1=max((currGC));
        currGC=1/scal1*currGC;
        
        %sort from left to right
        for i=1:length(Vxmin)
            if length(RES.peaks)>maxpeak_nr-1
                disp(['!! ABORT by reaching max peak number limit of: ' num2str(maxpeak_nr)])
                break
            end
            xmin=Vxmin(i);
            ymin=Vymin(i);
            xminidx=Vxminidx(i);
            
            % search max in close range
            j=xminidx-1;
            if j>=1
                while ymin<currGC(j)
                    ymin=currGC(j);
                    j=j-1;
                    if j<1
                        break
                    end
                end
            end
            if j~=xminidx-1
                xminidx=j+1;
                xmin=gx(xminidx);
            end
            
            j=xminidx+1;
            try
                if j<=length(gx)
                    while ymin<currGC(j)
                        ymin=currGC(j);
                        j=j+1;
                        if j>length(gx)
                            break
                        end
                    end
                end
            catch
                disp('<> Catch error: ID 0')
            end
            
            if j~=xminidx+1
                xminidx=j-1;
                xmin=gx(xminidx);
            end
            
            % approcimate peak width at 3/5 of max peak hight
            xidL=find(((gx<=xmin)+(   (currGC)<=ymin*3/5   ))==2,1,'last');
            widthL=xmin-gx(xidL);
            xidR=find(((gx>=xmin)+(   (currGC)<=ymin*3/5   ))==2,1,'first');
            widthR=gx(xidR)-xmin;
            abo=0;
            if isempty(widthL) || isempty(widthR)
                abo=1;
                widthid=1;
            else
                if widthL<widthR
                    xid=xidL;
                    width=widthL;
                else
                    xid=xidR;
                    width=widthR;
                end
                widthid=abs(xminidx-xid);
            end
            
 
            
            %ignore if peak too narrow
            if widthid<=2 || abo==1
                if info>1
                    disp(['<> continue (small peak): widthid=' num2str(widthid) '; abo=' num2str(abo)])
                end
                %ignore counter
                ign=ign+1;
                idxBlock(max(Vxminidx-ceil(para.bw*2),1):min(Vxminidx+floor(para.bw*2),length(idxBlock))  )=0;
                continue
            end
            
            % initial values
            p0=[0,...        center
                1.3*width,...       sigma left  (1.3*)
                1.3*width,...       sigma right (1.3*)
                0.95...          height
                0.3 ...           linearity left
                0.3 ...           linearity right
                ];
            
            % narrow peak range (fulfills: modell=profil)
            le_idx=max([1 xminidx-floor(0.5*widthid)]);
            ri_idx=min([length(gx) xminidx+ceil(0.5*widthid)]);
            
            % wide peak range (fulfills: modell<=profil)
            le_idx_w=max([1 xminidx-10*widthid]);
            ri_idx_w=min([length(gx) xminidx+10*widthid]);
            
            
            % Align center
            gxtemp=gx-xmin;
            actx=false(1,length(gx));
            actx(le_idx:ri_idx)=1;
            actx_w=false(1,length(gx));
            actx_w(le_idx_w:ri_idx_w)=1;
            
            % optimize
            lb=zeros(size(p0));
            ub=inf(size(p0));
            lb(1)=min(gxtemp(actx));
            ub(1)=max(gxtemp(actx));
            lb(2)=0.05*p0(2);
            lb(3)=0.05*p0(3);
            ub(2)=2.5*p0(2);
            ub(3)=10*p0(3);
            ub(5)=1;
            ub(6)=1;
            lb=lb(:);
            ub=ub(:);
            p0=p0(:);
            
            currtemp=currGC(actx);
            scaltemp=max(abs(currtemp));
            datf=(1/scaltemp)*currGC;
                      
            if lb(2)>=1E-5 && lb(3)>=1E-5
                if 1
                    we=[1 0];
                    fun2=@(p) tar_opt(p,gxtemp,datf,1,6,we,length(currtemp),length(currGC),actx,actx_w,0);
                    options=optimoptions(@lsqnonlin,...
                        'Display','off',...
                        'MaxIter',100,...
                        'MaxFunEvals',10000,...
                        'TolFun',1e-10,...
                        'TolX',1e-10);
                    pOpt=lsqnonlin(fun2,p0,lb,ub,options);

                else
                    pOpt=p0;
                end
                
                we=[1 10];
                fun2=@(p) tar_opt(p,gxtemp,datf,1,6,we,length(currtemp),length(currGC),actx,actx_w,0);
                options=optimoptions(@lsqnonlin,...
                    'Display','off',...
                    'MaxIter',100,...
                    'MaxFunEvals',10000,...
                    'FinDiffType','central',...
                    'TolFun',1e-10,...
                    'TolX',1e-10);
                pOpt=lsqnonlin(fun2,pOpt,lb,ub,options);
            else
                pOpt=p0;
            end
            
            if pOpt(2)<1E-8 || pOpt(3)<1E-8
                1;
            end
            
            
            % eval target function
            we=[1 1];
            fvec=tar_opt(pOpt,gxtemp,datf,1,6,we,length(currtemp),length(currGC),actx,actx_w,noiseLVL/scal1/scal/scaltemp);
            fev=sum(fvec.^2);
            
            pOpt(1)=pOpt(1)+xmin;
            pOpt(4)=scal1*scaltemp*pOpt(4);
            p0(1)=p0(1)+xmin;
            p0(4)=scal1*scaltemp*p0(4);
            
            
            mod=mod_ga(gx,pOpt);
                        
            % fit okay? (peak high enough, fit quality sufficient)
            if pOpt(4)>(noiseLVL/scal) && fev<fev_tol
                err=abs(trapz(gcs)-trapz(gcmod))/abs(trapz(gcs));
                if info>1
                    disp(['<> Add curve #' num2str(length(RES.peaks))])
                end
                
                if info>1
                    disp(['<> ' num2str(num) '; error: ' num2str(err) '; errortol: ' num2str(errTol)])
                end
                gcmod=gcmod+mod;
                
                st.mod=scal1;
                st.para=pOpt;
                RES.peaks{end+1}=st;
                curve_added=1;
            else
                if info>1
                    disp(['<> ignore fitted peak: hi= ' num2str(pOpt(4)) '<=' num2str((noiseLVL/scal)) '; fev= ' num2str(fev) '>=' num2str(fev_tol)])
                end
                ign=ign+1;
                idxBlock(max(Vxminidx-ceil(para.bw*2),1):min(Vxminidx+floor(para.bw*2),length(idxBlock))  )=0;
            end
        end
    else
        % If no peaks are detected, refine search parameters
        epsi = epsi/4;
        delta = delta/4;
    end
    if curve_added==0
        nocurve_prev=nocurve_prev+1;
    else
        nocurve_prev=0;
    end
end



%% Cluster Post Opt

P=[];
LB=[];
UB=[];

p_anz=0;
for i=1:length(RES.peaks)
    p0=RES.peaks{i}.para;
    if p0(1)>=gx(1) && p0(1)<=gx(end)
        p_anz=p_anz+1;
        P=[P; p0(:)];
        lb=zeros(size(p0));
        lb(2)=0.5*p0(2);
        lb(3)=0.5*p0(3);
        ub=inf(size(p0));
        ub(2)=2*p0(2);
        ub(3)=2*p0(3);
        ub(5)=1;
        ub(6)=1;
        LB=[LB; lb(:)];
        UB=[UB; ub(:)];
    end
end

Pm=reshape(P,6,p_anz);
LBm=reshape(LB,6,p_anz);
UBm=reshape(UB,6,p_anz);

%Aktive Bereiche
act=zeros(1,length(gx));
for i=1:p_anz
    spr=p_area;
    para=Pm(:,i);
    ce=para(1);
    siL=para(2);
    siR=para(3);
    act=act + (gx(:)'>=ce-spr*siL) - (gx(:)'>ce+spr*siR);
end
act=(act>=1);

targets=[];
area_ct=0;
for i=1:length(act)
    if i==1
        if act(i)==1
            area_ct=1;
            t1=1;
        end
    end
    if i>1
        if act(i-1)==0 && act(i)==1
            area_ct=area_ct+1;
            t1=i;
        end
        if act(i-1)==1 && act(i)==0
            targets=[targets;[t1 i]];
        end
    end
    if i==length(act)
        if act(i)==1
            targets=[targets;[t1 i]];
        end
    end
end

tar_zo=cell(area_ct,1);
for i=1:p_anz
    ce=Pm(1,i);
    for j=1:area_ct
        if  ce>=gx(targets(j,1)) && ce<=gx(targets(j,2))
            tar_zo{j}=[tar_zo{j} i];
            break
        end
        
    end
end



PP=zeros(size(Pm,1),length([tar_zo{:}]));
PR=cell(area_ct,1);

for i=1:area_ct
    pinit=Pm(:,tar_zo{i});
    Pt=pinit;
    PR{i}=reshape(Pt,6,length(tar_zo{i}));
end




for i=1:area_ct
    PP(:,tar_zo{i})=PR{i};
end

Pm=PP;

% update model
gcmod=zeros(length(gx),1);
peaks=cell(1,p_anz);
for i=1:p_anz
    peaks{i}.para=Pm(:,i);
    mod=mod_ga(gx,Pm(:,i));
    peaks{i}.mod=mod;
    gcmod=gcmod+mod;
end
RES.peaks=peaks;

%% Zuordnung kleine Kurven zu großen
% is a curve for more than p_corell*100 percent under another, it will be assigned to that one
if info>=2
    disp(' ')
    disp('=> Assignment')
end

Zuord=cell(length(RES.peaks),1);
for i=1:length(RES.peaks)
    Zuord{i}=i;
end


% overlap
for ii=1:overlap_runs
    if info>=2
        disp(['<> Overlap run ' num2str(ii) ' of ' num2str(overlap_runs)])
    end
    
    for i=1:length(RES.peaks)
        if ~isempty(Zuord{i})
            cZu=Zuord{i};
            mod=mod_ga(gx,Pm(:,cZu(1)));
            ce=Pm(1,cZu(1));
            siL=Pm(2,cZu(1));
            siR=Pm(3,cZu(1));
            leftb=ce-3*siL;
            rightb=ce+3*siR;
            for k=2:length(cZu)
                ce=Pm(1,cZu(k));
                siL=Pm(2,cZu(k));
                siR=Pm(3,cZu(k));
                leftb=min(leftb,  ce-3*siL );
                rightb=max(rightb,  ce+3*siR );
                mod=mod+mod_ga(gx,Pm(:,cZu(k)));
            end
        else
            continue
        end
        for j=1:length(RES.peaks)
            if i~=j
                if ~isempty(Zuord{j})
                    cZu2=Zuord{j};
                    
                    cce=Pm(1,cZu2(1));
                    
                    csiL=Pm(2,cZu2(1));
                    csiR=Pm(3,cZu2(1));
                    cleftb=cce-3*csiL;
                    crightb=cce+3*csiR;
                    for k=2:length(cZu2)
                        cce=Pm(1,cZu2(k));
                        csiL=Pm(2,cZu2(k));
                        csiR=Pm(3,cZu2(k));
                        cleftb=min(cleftb,  cce-3*csiL );
                        crightb=max(crightb,  cce+3*csiR );
                    end
                    
                    %fast check, trapz, interp1 etc cost time
                    if leftb>crightb || cleftb >rightb
                        continue
                    end
                    
                    mod2=mod_ga(gx,Pm(:,cZu2(1)));
                    for k=2:length(cZu2)
                        mod2=mod2+mod_ga(gx,Pm(:,cZu2(k)));
                    end
                    I_t=trapz(gx,mod);
                    I_t2=trapz(gx,mod2);
                    I_inter=trapz(gx,min(mod,mod2));
                    [m1,i1]=max(mod);
                    [m2,i2]=max(mod2);
                    modc=mod+mod2;
                    mc1=modc(i1);
                    mc2=modc(i2);
                    min_mc=min([mc1,mc2]);
                    mic=min(modc(min(i1,i2):max(i1,i2)));
 
                    % avoid peaks on tailing by skipping if valley between
                    % to maxima is deep enough
                    if abs((min_mc-mic)/min_mc)>=0.005 && ...
                            abs((I_t-I_inter)./I_inter)>0.005 && ...
                            abs((I_t2-I_inter)./I_inter)>0.005
                        continue
                    end
                    
                    if I_inter/I_t>p_correl || I_inter/I_t2>p_correl
                        Zuord{i}=[Zuord{i} Zuord{j}];
                        Zuord{j}=[];
                        continue
                    end
                    
                    
                    if i1<i2 && m1>m2*4 && I_inter/I_t2>0.1*p_correl && abs((mc2-mic)/mc2)<0.02
                        if info>1
                            disp('<> Tail found/appended')
                        end
                        Zuord{i}=[Zuord{i} Zuord{j}];
                        Zuord{j}=[];
                        continue
                    end
                    
                end
            end
        end
    end
end



RES.zuord=Zuord;
RES.gcmod=gcmod;














end

