classdef GCanalysis < hasHistory
    
    properties
        % raw data profile
        rawP profileC
        
        % baseline corrected data profile
        baseP profileC
        
        % clustered profile (peak regions devided by baseline regions)
        clusterP model
        
        % baseline after clustering, sparse
        basePS profileC
        
        % baseline after clustering, full
        basePF profileC
        
        % error measures
        modelError
        coveredPeakArea
        modelErrorNoN
        coveredPeakAreaNoN
        
        % progress indicator - contains time in seconds for the tasks
        progress=nan(6,1);
        
        % cluster result struct
        rescl
        
        % baseline result struct
        resbl
        
        % set of peaks
        peaks peak
        
        % signal after clustering, full
        profPF profileC
        
        % model after fit, full
        modePF profileC
        
        % algorithm settings
        para
        
        % timegrid range that is to be excluded
        % matrix(# x 2) - each row contains a left and right bound
        excludeRange=[]
        
        % timegrid range that is to be integrated as is
        % matrix(# x 2) - each row contains a left and right bound
        integrationRange=[]
        
        % peak-like class for range integration (length=size(integrationRange,1))
        intRangeArray=[]
        
        % system info - only relevant in validation mode
        sysInfo={};
        
        % parallel worker pool size to run threads on
        poolsize
        
        
    end
    
    properties(GetAccess=private)
        % labels for the values in GCanalysis.progress
        progresslabels={...
            'load data        ',...
            'initial baseline ',...
            'clustering       ',...
            'fitting          ',...
            'local corrections',...
            'total time       '}
    end
    
    
    methods
        %% Constructor
        function obj = GCanalysis(timegrid,signal,para)
            if nargin==2
                obj.rawP=profileC(timegrid,signal);
                lenref=ceil(gx(end)-gx(1));
                
                % default para
                para.info=1;                                    % Info-Level (0=off, 1=info, 2=debug, 3=trace)
                para.ver=-1;
                
                para.d=3;                                       % SG-Filter - poly degree
                para.M=10;                                      % SG-Filter - window half width 
                
                para.bl_n_all=max(ceil(1/10*lenref),10);        % Baseline - # support points for splines
                para.bl_relmaxpeakW=1/3;
                para.bl_startend=120;
                
                para.cl_n           =300;       % Cluster - # support points for clustering
                para.cl_first       =100;       % Cluster - # smallest value points for noise level approximation
                para.cl_fac         =1;         % Cluster - noise level scaling
                para.NoiseWin       =501;       % Cluster - MovingWindow width in noise level approx (has to be odd)
                
                para.maxpeak_nr     =min(ceil(4*lenref),350);    % Fit - # max peaks for each sub optimization problem
                para.fev_tol        =0.2;       % Fit - square sum tolerance in optimization
                para.errTol         =0.001;     % Fit - error tolerance tolerance in optimization
                para.retry          =50;
                para.epsi           =0.001;
                para.minPdist       =0.01;
                para.delta          =-1e-2;
                
                para.peakoverlap    =0.13;
                para.overlap_runs   =10;
                
                %Postprocessing
                para.relPeakHi      =1/6;
                para.relValley      =1/2;
                
                para.ispar          =1;
                
                para.bw = find((timegrid-timegrid(1))>para.minPdist,1,'first');
                
                para.isREF=0;
                para.dataR=[];
                
                para.tempf=[cd filesep 'GCana_temp' filesep];
                para.isknime=0;
                
                para.fName = '';
                
                obj.para=para;
                
            elseif nargin==3
                obj.rawP=profileC(timegrid,signal);
                obj.para=para;
            end
            
            mkdir(para.tempf)            
            obj.progress(1)=0;
            
        end
        
        % Destructor
        function delete(obj)
            p = properties(obj);
            for i=1:length(p)
                eval(['obj.' p{i} '=[];'])
            end
            clear obj
        end
        
        %% Auxillary
        %para
        function set.para(obj,para)
            obj.para=para;
        end
        
        function para=get.para(obj)
            para=obj.para;
        end
        
        %% Update peak integrals
        function updateIntegrals(obj)
            try
                % the fucntion adds information to each peak(i).intinfo
                % if already present, the information is updated
                
                % peak overlap range parameter
                nrDis=5;
                
                % check if peaks are present
                if isempty(obj.peaks)
                    warning('No peaks calculated. Call updateIntegrals after generatePeaks.')
                    return
                end
                
                % recalculate integral information
                gx=obj.rawP.timegrid;
                
                % get global signal area
                area_global=trapz(obj.profPF.timegrid,obj.profPF.signal);
                
                %sum of signal integrals
                sumSigAbs=0;
                sumAgiAbs=0;
                
                for i=1:length(obj.peaks)
                    %get current peak
                    cpeak=obj.peaks(i);
                    
                    %get model on full grid
                    fullpeak=cpeak.getmod(gx);
                    
                    %absolute model integral
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    intinfo.abs=trapz(gx,fullpeak);
                    
                    %relative model integral
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    intinfo.rel=intinfo.abs/area_global;
                    
                    
                    % get next nrDis peaks to the left and right in gx(left:right)
                    
                    %get peak bound
                    [ma,id]=max(fullpeak);
                    left=find(fullpeak>ma/300,1,'first');
                    right=find(fullpeak>ma/300,1,'last');
                    gx_l2r=gx(left:right);
                    
                    % get disjcunct peak idx range
                    disjunct_range = max(i-nrDis,1):min(i+nrDis,length(obj.peaks));
                    disjunct_range(disjunct_range==i)=[];
                    
                    % eval disjunct peaks
                    cMod_l2r=fullpeak(left:right);
                    disMod_l2r=zeros(length(gx_l2r),1);
                    for j=disjunct_range
                        disMod_l2r=disMod_l2r+obj.peaks(j).getmod(gx_l2r);
                    end
                    
                    % get ratio of current peak of all active peaks
                    area_cMod=trapz(gx_l2r,cMod_l2r);
                    area_disMod=trapz(gx_l2r,disMod_l2r);
                    ratio=area_cMod/(area_cMod+area_disMod);
                    
                    %absolute signal integral based on model bounds
                    intinfo.sigAbs=ratio*trapz(gx_l2r,obj.profPF.signal(left:right));
                    
                    %relative signal integral based on model bounds
                    % div by sum of signals is done subsequently
                    intinfo.sigRel=intinfo.sigAbs;
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    % Agilent integrals
                    
                    % Find global index of max of the peak
                    [~,max_id] = max(fullpeak);
                    
                    % Find global index of lbound of the peak
                    %l_id = find(fullpeak>=cpeak.signal(1),1,'first');
                    l_id = left;
                    
                    % Find global index of rbound of the peak
                    %r_id = find(fullpeak>=cpeak.signal(end),1,'last');
                    r_id = right;
                    
                    % Find global index of left minimum
                    [~,l_min_id] = min(fullpeak(l_id:max_id));
                    l_min_id = l_min_id+l_id-1;
                    
                    % Find global index of right minimum
                    [~,r_min_id] = min(fullpeak(max_id:r_id));
                    r_min_id = r_min_id+max_id-1;
                    
                    % Find the integral between the minimum
                    intinfo.agiAbs = trapz(gx(l_min_id:r_min_id),fullpeak(l_min_id:r_min_id));
                    
                    %relative signal integral based on model bounds
                    % div by sum of signals is done subsequently
                    intinfo.agiRel=intinfo.agiAbs;
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % write to peak
                    obj.peaks(i).intinfo=intinfo;
                    
                    sumSigAbs=sumSigAbs+intinfo.sigAbs;
                    sumAgiAbs=sumAgiAbs+intinfo.agiAbs;
                end
                
                for i=1:length(obj.peaks)
                    obj.peaks(i).intinfo.sigRel=obj.peaks(i).intinfo.sigRel/sumSigAbs;
                    obj.peaks(i).intinfo.agiRel=obj.peaks(i).intinfo.agiRel/sumAgiAbs;
                end
            catch ex
                obj.errorLog(ex,'update integrals')
            end
        end
        
        %% Get peak integrals
        function [inteAbsMod, inteRelMod, inteAbsSig, inteRelSig, inteAbsAgi, inteRelAgi] = getIntegrals(obj)
            try
                % check if peaks are present
                if isempty(obj.peaks)
                    warning('No peaks calculated. Call updateIntegrals after generatePeaks.')
                    return
                end
                
                %absolute int model based
                inteAbsMod = zeros(length(obj.peaks),1);
                %relative int model based
                inteRelMod = zeros(length(obj.peaks),1);
                
                %absolute int signal based
                inteAbsSig = zeros(length(obj.peaks),1);
                %relative int signal based
                inteRelSig = zeros(length(obj.peaks),1);
                
                %absolute int agilent
                inteAbsAgi = zeros(length(obj.peaks),1);
                %relative int agilent
                inteRelAgi = zeros(length(obj.peaks),1);
                
                for i=1:length(obj.peaks)
                    inteAbsMod(i)=obj.peaks(i).intinfo.abs;
                    inteRelMod(i)=obj.peaks(i).intinfo.rel;
                    inteAbsSig(i)=obj.peaks(i).intinfo.sigAbs;
                    inteRelSig(i)=obj.peaks(i).intinfo.sigRel;
                    inteAbsAgi(i)=obj.peaks(i).intinfo.agiAbs;
                    inteRelAgi(i)=obj.peaks(i).intinfo.agiRel;
                end
                
            catch ex
                obj.errorLog(ex,'get peak integrals')
            end
        end
        
        
        %% Calculation Baseline
        function calcBaseline(obj,meth,BLopt)
            try
                if obj.para.info>=1
                    disp(['### Baseline calc - Start ###'])
                end
                
                % BL complete time axis
                tic
                if nargin==1
                    meth='new';
                end
                if ischar(meth)
                    if strcmp(meth,'new') || strcmp(meth,'beads')
                    else
                        meth='new';
                        warning('Baseline method unknown: continue with standard')
                    end
                else
                    meth='new';
                    warning('Baseline method unknown: continue with standard')
                end
                
                switch meth
                    case 'new'
                        z=obj.getBaseline();
                    case 'beads'
                        z=obj.getBaselineBEADS(BLopt);
                end
                
                sig=obj.rawP.signal;
                
                obj.progress(2)=toc;
                obj.progress(6)=sum(obj.progress(1:5));
                obj.baseP=profileC(obj.rawP.timegrid,sig-z);
                
                if obj.para.info==2
                    save([obj.para.tempf 'baseline'])
                end
                if obj.para.info>=1
                    disp(['### Baseline calc - End ###'])
                end
            catch ex
                obj.errorLog(ex,'calculate baseline')
            end
        end
        
        
        function z=getBaselineBEADS(obj,BLopt)
            try
                sig=obj.rawP.signal;
                timeg=obj.rawP.timegrid;
                z=zeros(size(sig));
            catch ex
                obj.errorLog(ex,'getBaselineBEADS')
            end
        end
        
        
        
        function z=getBaseline(obj)
            try
                % maximale peakbreite in einem Unterteilungsschritt
                % wenn der breiteste peak mehr als relmaxpeak des
                % GC-Ausschnitts einnimmt wird nicht weiter unterteilt.
                relmaxpeakW=obj.para.bl_relmaxpeakW;
                PP_relmaxpeakW=obj.para.bl_PP_relmaxpeakW;
                
                %lade raw gc
                sig=obj.rawP.signal;
                timeg=obj.rawP.timegrid;
                
                % initialisiere Unterteilung
                % aqGr - Indizes der Unterteilung
                aqGr=[1 length(timeg)];
                % aqAc - Aktive noch zu analysierende Unterteilungsabschnitte
                % aqAc(i) -> Bereich [i i+1]
                aqAc=1;
                while max(aqAc)==1
                    % aqGrAdd - Hinzuzufuegende Indizes der Unterteilung
                    aqGrAdd=[];
                    % finde aktive Abschnitte
                    aqid=find(aqAc==1);
                    % init alle Abschnitte als inaktiv
                    aqAc=zeros(length(aqGr)-1,1);
                    
                    % pruefe alle vormals aktiven Abschnitte
                    for i=1:length(aqid)
                        % Unterteilungsabschnitt bestimmen
                        % get signal
                        aqSig=sig(aqGr(aqid(i)):aqGr(aqid(i)+1));
                        % get grid
                        aqTime=timeg(aqGr(aqid(i)):aqGr(aqid(i)+1));
                        
                        % Unterteilungsabschnitt weiter unterteilen?
                        if isdevide(aqTime,aqSig,100,relmaxpeakW)
                            %mittlerer Index
                            midtimeidx=round(length(aqTime)/2);
                            %zugehoerige Zeit
                            mid=aqTime(midtimeidx);
                            
                            % get total, start and end time + indices
                            % the section covers a range around the middle of
                            % the time grid, where the divide index is
                            % searched.
                            % take just the direct middle can "hit" a peak and
                            % results in wrong baselines
                            timelen=aqTime(end)-aqTime(1);
                            starttime=(mid-0.5*relmaxpeakW*timelen);
                            endtime=(mid+0.5*relmaxpeakW*timelen);
                            startidx=find(aqTime>=starttime,1,'first');
                            endidx=find(aqTime<=endtime,1,'last');
                            % section to small?
                            if endidx<=startidx+20
                                break
                            end
                            
                            % linear baseline substraction based on signal start
                            % and end
                            Cor=linspace(aqSig(1),aqSig(end),length(aqTime));
                            aqSigLinCor=aqSig(:)-Cor(:);
                            
                            % get minimum in inner range
                            [~,testid]=min(aqSigLinCor(startidx:endidx));
                            % devide index
                            minidx=startidx-1+testid;
                            
                            aqGrAdd(end+1)=aqGr(aqid(i))-1+minidx;
                            
                            aqAc(aqid(i))=1;
                        end
                    end
                    aqGrUnsort=[aqGr(1:end-1); aqGrAdd(:)];
                    aqAc=[aqAc(:); ones(length(aqGrAdd),1)];
                    
                    [aqGr,id]=sort(aqGrUnsort);
                    aqGr(end+1)=length(timeg);
                    aqAc=aqAc(id);
                end
                
                % Postproc - negative Bereiche eleminieren
                % more baseline points are added/points are adapted, where a
                % negative corrected signal would occur.
                
                
                aqGrPre=aqGr;
                
                %get baselone approximation
                [intx,inty]=grid2val(obj,aqGr);
                if length(intx)==1
                    z=min(sig)*ones(size(timeg));
                else
                    z=pchip(intx,inty,timeg);
                end
                
                %nr post process runs
                bl_PPruns=obj.para.bl_PPruns;
                for r=1:bl_PPruns
                    %nr inserted bl points
                    nrinsterted=0;
                    aqGrAdd=[];
                    for i=1:length(aqGr)-1
                        cTime=timeg(aqGr(i):aqGr(i+1)-1);
                        cSig=sig(aqGr(i):aqGr(i+1)-1);
                        cBL=z(aqGr(i):aqGr(i+1)-1);
                        cCor=cSig-cBL;
                        sigIntTot=abs(trapz(cTime,cCor));
                        sigIntNeg=abs(trapz(cTime,min(cCor,0)));
                        
                        % check ration of negative area against total area
                        % if true, add another baseline point analogous to the
                        % addition of a divide index above.
                        if sigIntNeg/sigIntTot>0.1
                            nrinsterted=nrinsterted+1;
                            aqSig=cSig;
                            aqTime=cTime;
                            
                            midtimeidx=round(length(aqTime)/2);
                            mid=aqTime(midtimeidx);
                            
                            timelen=aqTime(end)-aqTime(1);
                            starttime=(mid-0.5*PP_relmaxpeakW*timelen);
                            endtime=(mid+0.5*PP_relmaxpeakW*timelen);
                            startidx=find(aqTime>=starttime,1,'first');
                            endidx=find(aqTime<=endtime,1,'last');
                            if endidx<=startidx+20
                                break
                            end
                            [~,testid]=min(aqSig(startidx:endidx));
                            minidx=startidx-1+testid;
                            
                            addI=aqGr(i)-1+minidx;
                            aqGrAdd(end+1)=addI;
                            
                            if 0
                                figure(34);
                                clf
                                plot(cTime,cSig)
                                hold on
                                plot(cTime,cBL)
                                plot(cTime(minidx),0,'r*')
                            end
                        end
                    end
                    aqGr=sort([aqGr(:); aqGrAdd(:)]);
                    if obj.para.info>=2
                        disp(['=> Run ' num2str(r) ', Inserted Points: ' num2str(nrinsterted)])
                    end
                    
                    %update BL
                    [intx,inty]=grid2val(obj,aqGr);
                    if length(intx)==1
                        z=min(sig)*ones(size(timeg));
                    else
                        z=pchip(intx,inty,timeg);
                    end
                end
                
                
                if obj.para.info>=3
                    figure
                    plot(intx,inty,'r*')
                    hold on
                    plot(timeg,z,'r')
                    plot(obj.rawP,'k')
                    plot(obj.rawP.timegrid(aqGrPre),zeros(size(aqGrPre)),'g*')
                    legend('selected data','bl','prof raw','interval borders')
                    
                end
            catch ex
                obj.errorLog(ex,'getBaselineDef')
            end
        end
        
        %% calc bl from grid points
        function [intx,inty]=grid2val(obj,aqGr)
            % Hinzufuegen des ersten und letzten Punktes im GC zur BL.
            % y-Wert wird auf Basis des Minimums der ersten oder letzten
            % "startend" Punkte bestimmt.
            try
                startend=obj.para.bl_startend;
                
                sig=obj.rawP.signal;
                timeg=obj.rawP.timegrid;
                localmmean=zeros(length(aqGr)-1,1);
                
                C_act_x2=zeros(length(aqGr)-1,1);
                for i=1:length(aqGr)-1
                    [C_D2,miIdx]=min(sig(aqGr(i):aqGr(i+1)-1));
                    localmmean(i)=C_D2;
                    
                    C_act_x2(i)=timeg(aqGr(i)-1+miIdx);
                end
                
                intx=C_act_x2(:);
                inty=localmmean(:);
                
                if intx(1)~=timeg(1)
                    intx=[timeg(1); intx];
                    inty=[min(sig(1:startend+1)); inty];
                end
                
                if intx(end)~=timeg(end)
                    intx=[intx; timeg(end)];
                    inty=[inty; min(sig(end-startend:end))];
                end
                
            catch ex
                obj.errorLog(ex,'grid2val')
            end
        end
        
        %% Cluster bestimmen
        function calcCluster(obj)
            try
                if obj.para.info>=1
                    disp(['### Cluster calc - Start ###'])
                end
                tic
                %apply cluster algorithm (see bin/corr_cluster.m)
                obj.rescl=obj.baseP.calcClusterP(obj.para);
                %create individual cluster profiles
                obj.clusterP(obj.rescl.anz_prob,1)=profileC();
                
                switch 'no BL'
                    %                 case 'local BL'
                    %                 info=obj.para.info;
                    %
                    %                     temp_clusterP{obj.rescl.anz_prob,1}=model();
                    %                     for i=1:obj.rescl.anz_prob
                    %                         temp_clusterP{i}=model(obj.rescl.Cgx{i},obj.rescl.Cgc{i});
                    %                     end
                    %
                    %                     resbl_temp=cell(obj.rescl.anz_prob,1);
                    %                     para_temp=obj.para;
                    %
                    %                     %%%
                    %                     if obj.para.ispar==1
                    %                         parfor i=1:obj.rescl.anz_prob
                    %                             resbl_temp{i}=temp_clusterP{i}.calcBaselineP(para_temp,8);
                    %                             if info>=2
                    %                                 disp(['=> Finished baseline cluster ' num2str(i)]);
                    %                             end
                    %                         end
                    %                     else
                    %                         for i=1:obj.rescl.anz_prob
                    %                             resbl_temp{i}=temp_clusterP{i}.calcBaselineP(para_temp,8);
                    %                             if info>=2
                    %                                 disp(['=> Finished baseline cluster ' num2str(i)]);
                    %                             end
                    %                         end
                    %                     end
                    
                    case 'no BL'
                        resbl_temp=cell(obj.rescl.anz_prob,1);
                        for i=1:obj.rescl.anz_prob
                            resbl_temp{i}.bl.mod=zeros(size(obj.rescl.Cgx{i}));
                            resbl_temp{i}.gcs=obj.rescl.Cgc{i};
                        end
                end
                
                % combine individual cluster information to full grid info
                bl=[];
                gxr=[];
                for i=1:obj.rescl.anz_prob
                    obj.clusterP(i)=model(obj.rescl.Cgx{i},resbl_temp{i}.gcs);
                    bl=[bl; resbl_temp{i}.bl.mod];
                    gxr=[gxr; obj.rescl.Cgx{i}];
                end
                %baseline, sparse
                obj.basePS=profileC(gxr,bl);
                %baseline, full grid
                obj.basePF=profileC(obj.rawP.timegrid,interp1(gxr,bl,obj.rawP.timegrid,'linear',0));
                %corrected signal, full grid
                obj.profPF=profileC(obj.rawP.timegrid,obj.baseP.signal-obj.basePF.signal);
                
                obj.progress(3)=toc;
                obj.progress(6)=sum(obj.progress(1:5));
                
                if obj.para.info==2
                    save([obj.para.tempf 'cluster'])
                end
                if obj.para.info>=1
                    disp(['### Cluster calc - End ###'])
                end
            catch ex
                obj.errorLog(ex,'calcCluster')
            end
        end
        
        %% Start parallel matlab workers pool
        function startParPool(obj)
            try
                nW=obj.para.numWorkers;
%                 if nW==1
%                     obj.para.ispar=0;
%                     warning('parallel computing pool has not been started. Too less cpu cores. [>=3 are recommended]')
%                 end
                
                %check if pool is open
                p = gcp('nocreate');
                if ~isempty(p)
                    if p.NumWorkers==nW
                        warning('ParPool already open with correct number of workers')
                        obj.poolsize=nW;
                        return
                    else
                        delete(p);
                    end
                end
                
                if obj.para.ispar==1
                    try
                        parpool(nW)
                    catch
                        error('gcp could not be started.')
                    end
                else
                    if ~isempty(gcp('nocreate'))
                        delete(gcp)
                    end
                end
                
                % verify active number of workers
                p = gcp('nocreate');
                if isempty(p)
                    obj.para.ispar=0;
                    ps = 0;
                else
                    ps = p.NumWorkers;
                end
                obj.poolsize=ps;
                
            catch ex
                obj.errorLog(ex,'startParPool')
            end
        end
        
        %% End parallel matlab workers pool
        function endParPool(obj)
            if ~isempty(gcp('nocreate'))
                delete(gcp)
            end
        end
        
        
        %% Fit Modell Parallel (HAUPTTEIL)
        function calcFitPar(obj)
            try
                if obj.para.info>=1
                    disp(['### Fit calc par - Start ###'])
                end
                tic
                
                
                %get noise level
                noiselvl=obj.rescl.noiselvl;
                
                %maximum number of queued jobs above poolsize
                pooloverload=obj.para.pooloverload;
                
                ps=obj.poolsize;
                jobs(1:ps+pooloverload) = parallel.FevalFuture;
                
                % indicates if all computations in a cluster are finished
                clusterAct=ones(size(obj.clusterP));
                %indicates zones on the time grid in each cluster for which a
                %calculation is runnung (0 = no calc, 1 = calc)
                clusterCurCalcZone=cell(size(obj.clusterP));
                
                for i=1:length(clusterCurCalcZone)
                    clusterCurCalcZone{i}=zeros(size(obj.clusterP(i).timegrid));
                end
                
                % next cluster to retrieve a peak calc from
                ClusterIdx=1;
                
                %currently used workers
                %             currps=0;
                
                % counter for failed peak init determination
                NextPeakCounter=0;
                
                % jobinfo - i-th entry indicates, if a job is assigned
                %           to the i-th matlab worker
                jobinfo=zeros(1,ps+pooloverload);
                
                %cluster assignment - i-th entry with value x indicates:
                % x=0  -> no current job on matlab worker
                % x>0  -> i-th matlab worker runs a job for cluster x
                clusterinfo=zeros(1,ps+pooloverload);
                
                isDetailedPlot=obj.para.isDetailedPlot;%obj.rescl.anz_prob;
                %update details plot every "updaterate" loop iterations
                updaterate=50;
                if isDetailedPlot
                    f=figure(99);
                    
                    % fit error
                    PlErr=zeros(length(obj.clusterP),1);
                    for i=1:length(obj.clusterP)
                        res=obj.clusterP(i).signal-obj.clusterP(i).modeleval;
                        err=abs(trapz(res))/abs(trapz(obj.clusterP(i).signal));
                        PlErr(i,1)=err;
                    end
                    
                    % #peak
                    PlNr=zeros(length(obj.clusterP),1);
                    for i=1:length(obj.clusterP)
                        PlNr(i,1)=length(obj.clusterP(i).peaks);
                    end
                end
                
                % outer loop (active if clusters are active or jobs running)
                ctn=1;
                while max(clusterAct)==1 || max(jobinfo)==1
                    
                    % update info
                    %%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if isDetailedPlot
                        ctn=ctn+1;
                        if updaterate==ctn
                            ctn=1;
                            
                            % append current fit error
                            addv=zeros(length(obj.clusterP),1);
                            for i=1:length(obj.clusterP)
                                res=obj.clusterP(i).signal-obj.clusterP(i).modeleval;
                                err=abs(trapz(res))/abs(trapz(obj.clusterP(i).signal));
                                addv(i,1)=err;
                            end
                            PlErr=[PlErr addv];
                            
                            % append current peak number
                            for i=1:length(obj.clusterP)
                                addv(i,1)=length(obj.clusterP(i).peaks);
                            end
                            PlNr=[PlNr addv];
                            
                            % obtain still active clusters
                            lc=cell(length(obj.clusterP),1);
                            for i=1:length(obj.clusterP)
                                if clusterAct(i)
                                    lc{i}='active';
                                else
                                    lc{i}='off';
                                end
                            end
                            
                            %plot
                            figure(f)
                            clf
                            subplot(1,2,1)
                            semilogy(PlErr');
                            legend(lc,'Location','SouthWest')
                            title(['err abs(trapz(sig-mod))/abs(trapz(sig)), limit=' num2str(obj.para.errTol) ])
                            
                            subplot(1,2,2)
                            plot(PlNr');
                            title(['Added peaks, max=' num2str(obj.para.maxpeak_nr)])
                            drawnow
                        end
                    end
                    
                    % check for finished jobs
                    %%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    for i=1:ps+pooloverload
                        
                        % get job state
                        state=jobs(i).State;
                        
                        % check for failed jobs, if present -> error
                        if strcmp(state,'failed')
                            save test_failed
                            error('job failed')
                        end
                        
                        % check for finished unread jobs
                        if strcmp(state,'finished') && ~jobs(i).Read
                            if obj.para.info>=2
                                disp(['par <> job ' num2str(jobs(i).ID) ' finished on idx ' num2str(i)])
                            end
                            
                            %update monitor variables
                            jobinfo(i)=0;
                            clusterinfo(i)=0;
                            
                            % get results
                            resJob=jobs(i).fetchOutputs;
                            
                            %cluster idx in clusterP
                            outCluIdx=resJob.ClusterIdx;
                            
                            %update calcWin
                            clusterCurCalcZone{outCluIdx}=clusterCurCalcZone{outCluIdx}-resJob.calcWin;
                            
                            init_idx=resJob.init_idx; %timegridIdx in job creation
                            
                            %check for negativity
                            mod=resJob.mod; %unscaled model on full cluster time grid
                            newmod=obj.clusterP(outCluIdx).modeleval+mod;
                            sig=obj.clusterP(outCluIdx).signal;
                            if min(sig(resJob.calcWin)-newmod(resJob.calcWin))<=-0.1*max(sig(resJob.calcWin))
                                if obj.para.info>=1
                                    disp(['par <> peak fit failed (negativity)'])
                                end
                                resJob.success=0;
                            end
                            
                            %check if res peak is in exlcude zone
                            for i_ex=1:size(obj.excludeRange,1)
                                if resJob.pOpt(1)>=obj.excludeRange(i_ex,1) && resJob.pOpt(1)<=obj.excludeRange(i_ex,2)
                                    resJob.success=0;
                                    % add history
                                    if obj.isHistory
                                        hist_comment=[];
                                        hist_comment.jobidx = i;
                                        hist_comment.deadzone = obj.clusterP(outCluIdx).deadzone;
                                        hist_comment.success = resJob.success;
                                        hist_comment.cluster = outCluIdx;
                                        hist_comment.calcWindow = resJob.calcWin;
                                        hist_comment.clusterAct = clusterAct(outCluIdx);
                                        if resJob.success
                                            hist_comment.pOpt = resJob.pOpt;
                                        end
                                        
                                        obj.addHist(hist_comment)
                                    end
                                    break
                                end
                            end
                            
                            
                            if resJob.success==0
                                if obj.para.info>=2
                                    disp(['par <> peak fit failed (fit error)'])
                                end
                                %if peak fit fails, the cluster deadzone is extended
                                obj.clusterP(outCluIdx).deadzone(max(init_idx-ceil(obj.para.bw*2),1):min(init_idx+floor(obj.para.bw*2),length(obj.clusterP(outCluIdx).signal))  )=false;
                                obj.clusterP(outCluIdx).noCurveAddCounter=obj.clusterP(outCluIdx).noCurveAddCounter+1;
                                % check if cluster stays active (if not already inactive)
                                if clusterAct(outCluIdx)
                                    clusterAct(outCluIdx)=obj.clusterP(outCluIdx).checkClusterAct(obj.para,outCluIdx,noiselvl);
                                end
                                % add history
                                if obj.isHistory
                                    hist_comment=[];
                                    hist_comment.jobidx = i;
                                    hist_comment.deadzone = obj.clusterP(outCluIdx).deadzone;
                                    hist_comment.success = resJob.success;
                                    hist_comment.cluster = outCluIdx;
                                    hist_comment.calcWindow = resJob.calcWin;
                                    hist_comment.clusterAct = clusterAct(outCluIdx);
                                    if resJob.success
                                        hist_comment.pOpt = resJob.pOpt;
                                    end
                                    
                                    obj.addHist(hist_comment)
                                end
                                
                                % => retry at another spot in the same cluster
                                continue
                            end
                            
                            % counts how many fits have failed in a row (cf. "% check cluster activity" below)
                            obj.clusterP(outCluIdx).noCurveAddCounter=0;
                            
                            % postprocessing for successful peak fit
                            
                            % update model in cluster
                            obj.clusterP(outCluIdx).modeleval=newmod;
                            
                            % add peak
                            sp.para=resJob.pOpt;
                            sp.mod=mod;
                            obj.clusterP(outCluIdx).peaks{end+1}=sp;
                            
                            % check if cluster stays active (if not already inactive)
                            if clusterAct(outCluIdx)
                                clusterAct(outCluIdx)=obj.clusterP(outCluIdx).checkClusterAct(obj.para,outCluIdx,noiselvl);
                            end
                            
                            % add history
                            if obj.isHistory
                                hist_comment=[];
                                hist_comment.jobidx = i;
                                hist_comment.deadzone = obj.clusterP(outCluIdx).deadzone;
                                hist_comment.success = resJob.success;
                                hist_comment.cluster = outCluIdx;
                                hist_comment.calcWindow = resJob.calcWin;
                                hist_comment.clusterAct = clusterAct(outCluIdx);
                                if resJob.success
                                    hist_comment.pOpt = resJob.pOpt;
                                end
                                
                                obj.addHist(hist_comment)
                            end
                            
                        end
                    end
                    if min(jobinfo)==1
                        %idle if all workers are busy
                        pause(0.02)
                        %check again
                        continue
                    end
                    
                    % add jobs
                    %%%%%%%%%%
                    breakw=0;
                    
                    while sum(jobinfo)<ps+pooloverload
                        if breakw
                            break
                        end
                        breakw=0;
                        
                        %check if active clusters exist
                        if max(clusterAct)==0
                            break
                        end
                        %skip inactive cluster
                        while clusterAct(ClusterIdx)==0
                            ClusterIdx=ClusterIdx+1;
                            if ClusterIdx > length(clusterAct)
                                ClusterIdx=1;
                            end
                        end
                        
                        % check for max calculations per cluster
                        currCluster=obj.clusterP(ClusterIdx);
                        if sum(clusterinfo==ClusterIdx) < obj.para.maxCalcPerClus
                            % get next peak to calc outside of deadzone and currCalcZone
                            [timegridIdx,smooth_sig_max]=currCluster.getNextPeakCalc(clusterCurCalcZone{ClusterIdx},obj.para);
                        else
                            % maximum number of calculations per cluster
                            % reached, continue with next cluster
                            timegridIdx=-1;
                        end
                        
                        if timegridIdx==-1
                            % check if cluster stays active (if not already inactive)
                            if clusterAct(ClusterIdx)
                                clusterAct(ClusterIdx)=currCluster.checkClusterAct(obj.para,ClusterIdx,noiselvl);
                            end
                            
                            NextPeakCounter=NextPeakCounter+1;
                            if sum(clusterAct)==NextPeakCounter
                                % pause if no peak position found in all
                                % remaining free clusters
                                breakw=1;
                                NextPeakCounter=0;
                            else
                                %otherwise try next cluster
                                ClusterIdx=ClusterIdx+1;
                                if ClusterIdx > length(clusterAct)
                                    ClusterIdx=1;
                                end
                            end
                            continue
                        else
                            NextPeakCounter=0;
                        end
                        
                        % get approximate peak infos
                        %p => struct with fields center, width, range, etc
                        p=currCluster.getApproxPeak(timegridIdx,smooth_sig_max,obj.para);
                        if p.success==0
                            %if peak initialzation fails, the cluster deadzone
                            %is extended (inside getApproxPeak)
                            % => retry at another spot in the same cluster
                            continue
                        end
                        
                        % generate a peak fit problem
                        j=currCluster.prepareJob(p,noiselvl);
                        
                        %update calcWin
                        clusterCurCalcZone{ClusterIdx}=clusterCurCalcZone{ClusterIdx}+j.calcZone(:);
                        
                        %get free job idx
                        jobidx=find(jobinfo==0,1,'first');
                        
                        % send prob to pool / solve
                        jobs(jobidx)=parfeval(@job_fcn,1,j,p,currCluster.timegrid,noiselvl,ClusterIdx,obj.para.fev_tol,obj.para.info,timegridIdx,j.calcZone(:));
                        
                        
                        jobinfo(jobidx)=1;
                        clusterinfo(jobidx)=ClusterIdx;
                        
                        if obj.para.info>=2
                            try
                                disp(['par <> job ' num2str(jobs(jobidx).ID) ' pushed'])
                            catch
                                save temp_err
                                error('aa')
                            end
                        end
                        
                        % add history
                        if obj.isHistory
                            hist_comment=[];
                            hist_comment.jobidx = jobidx;
                            hist_comment.deadzone = obj.clusterP(ClusterIdx).deadzone;
                            hist_comment.success = -1; % job not started
                            hist_comment.cluster = ClusterIdx;
                            hist_comment.calcWindow = clusterCurCalcZone{ClusterIdx};
                            hist_comment.clusterAct = clusterAct(ClusterIdx);
                            hist_comment.p = p;
                            
                            obj.addHist(hist_comment)
                        end
                        
                        ClusterIdx=ClusterIdx+1;
                        if ClusterIdx > length(clusterAct)
                            ClusterIdx=1;
                        end
                    end
                end
                
                %convert to old format
                for i=1:length(obj.clusterP)
                    obj.clusterP(i).resfi.peaks=obj.clusterP(i).peaks;
                    obj.clusterP(i).resfi.scal=1;
                    obj.clusterP(i).resfi.zuord=obj.clusterP(i).calcZuordAlt2(obj.para);
                end
                
                
                %init/update model
                gxm=[];
                for i=1:length(obj.clusterP)
                    gxm=[gxm; obj.clusterP(i).modeleval];
                end
                obj.modePF=profileC(obj.rawP.timegrid, interp1(obj.basePS.timegrid ,gxm,obj.basePF.timegrid,'linear',0));
                
                
                obj.progress(4)=toc;
                obj.progress(6)=sum(obj.progress(1:5));
                if obj.para.info==2
                    save([obj.para.tempf 'fit'])
                end
                if obj.para.info>=1
                    disp(['### Fit calc par - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'calcFit')
            end
        end
        
        
        %% Correct local high error
        function obj=correctLocalError(obj,runs)
            try
                if obj.para.info>=1
                    disp(['### Local error correction - Start ###'])
                end
                tic
                gx=obj.profPF.timegrid;
                ma_global=max(obj.profPF.signal);
                
                E=obj.getError;
                
                E0=norm(E,2);
                for rr=1:runs
                    if obj.para.info>=2
                        disp(['=> Run ' num2str(rr) ' of ' num2str(runs)])
                    end
                    
                    [~,idxEmax]=max(E);
                    startidx=obj.getClosestPeak2IDX(idxEmax);
                    
                    I=obj.getConPeaks(0.01,startidx);
                    
                    
                    % get bounds
                    lb=obj.peaks(startidx).peakstruct.lbound;
                    rb=obj.peaks(startidx).peakstruct.rbound;
                    for i=I
     
                        lbn=obj.peaks(i).peakstruct.lbound;
                        rbn=obj.peaks(i).peakstruct.rbound;
                        if lbn<=lb
                            lb=lbn;
                        end
                        if rbn>=rb
                            rb=rbn;
                        end
                    end
                    idxL=find(gx>=lb,1,'first');
                    idxR=find(gx<=rb,1,'last');
                    
                    E(idxL:idxR)=0;
                    obj=obj.postLocalFit(startidx);
                    
                    for i=I
                        obj.peaks(i).updatePropAfterCorr(gx,ma_global);
                    end
                end
                
                obj.updateModelAndError;
                
                if obj.para.info>=2
                    E=obj.getError;
                    Ee=norm(E,2);
                    disp(['=> Results'])
                    disp(['   Error before: ' num2str(E0) '; Error after: ' num2str(Ee)])
                end
                
                obj.progress(5)=toc;
                obj.progress(6)=sum(obj.progress(1:5));
                
                if obj.para.info>=1
                    disp(['### Local error correction - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'correctLocalError')
            end
        end
        
        %% Update model and model error
        function obj=updateModelAndError(obj)
            try
                gx=obj.profPF.timegrid;
                mod=zeros(size(gx));
                for i=1:length(obj.peaks)
                    mod=mod+obj.peaks(i).getmod(gx);
                end
                obj.modePF.signal=mod;
                
                
                %remove exclude ranges (no arg=full range)
                exIdxRange=obj.getExcludeIndices;
                sig=obj.profPF.signal(:);
                sig=sig(exIdxRange);
                mod=mod(exIdxRange);
                
                
                %Update errors
                obj.modelError=norm(sig-mod,2)/norm(sig(:),2);
                obj.coveredPeakArea=1-trapz(abs(sig(:)-mod(:)))/trapz(abs(sig(:)));
                
                noiselvl=obj.rescl.noiselvl;
                obj.modelErrorNoN=norm(max(abs(sig(:)-mod(:))-noiselvl,0)  ,2)/norm(sig(:),2);
                obj.coveredPeakAreaNoN=1-trapz(max(abs(sig(:)-mod(:))-noiselvl,0))/trapz(abs(sig(:)));
                
            catch ex
                obj.errorLog(ex,'updateModelAndError')
            end
        end
        
        %% Generate Peaks from calc result
        function obj=generatePeaks(obj)
            try
                if obj.para.info>=1
                    disp(['### Generate peaks - Start ###'])
                end
                ma_global=max(obj.profPF.signal);
                fix=1;
                PEAKS={};
                for jj=1:obj.rescl.anz_prob
                    if obj.para.info>=2
                        disp(['=> Cluster # ' num2str(jj)]);
                    end
                    allRES=obj.clusterP(jj).resfi;
                    scal=allRES.scal;
                    currX=obj.rescl.Cgx{jj};                    
                    
                    
                    for i=1:length(allRES.peaks)
                        zuord=allRES.zuord{i};
                        
                        if ~isempty(zuord)
                            paramat=zeros(length(allRES.peaks{zuord(1)}.para),length(zuord));
                            
                            if obj.para.info>=3
                                disp(['   Peak # ' num2str(i) ', ' num2str(fix)]);
                            end
                            fullpeak=allRES.peaks{zuord(1)}.mod*scal;
                            
                            paramat(:,1)=allRES.peaks{zuord(1)}.para(:);
                            
                            for j=2:length(zuord)
                                fullpeak=fullpeak+allRES.peaks{zuord(j)}.mod*scal;
                                paramat(:,j)=allRES.peaks{zuord(j)}.para(:);
                            end
                            paramat(4,:)=scal*paramat(4,:);
                            
                            %sort paramat
                            [~,idsort]=sort(paramat(1,:));
                            paramat=paramat(:,idsort);
                            
                            
                            [ma,id]=max(fullpeak);
                            left=find(fullpeak>ma/300,1,'first');
                            right=find(fullpeak>ma/300,1,'last');
                            if isempty(left)
                                left=1;
                            end
                            if isempty(right)
                                right=length(fullpeak);
                            end
                            inte=trapz(currX,fullpeak);
                            if obj.para.info>=3
                                disp(['     Left bound:        ' num2str(currX(left))]);
                                disp(['     Peak max location: ' num2str(currX(id))]);
                                disp(['     Right bound:       ' num2str(currX(right))]);
                                disp(['     Peak height:       ' num2str(ma)]);
                                disp(['     Peak rel. height:  ' num2str(ma/ma_global)]);
                                disp(['     Model number:      ' num2str(length(zuord))]);
                            end
                            p.peakX_l2r=currX(left:right);
                            p.fullpeak_l2r=fullpeak(left:right);
                            p.paramat=paramat;
                            p.lbound=currX(left);
                            p.rbound=currX(right);
                            p.maxlocation=currX(id);
                            p.height=ma;
                            p.relheight=ma/ma_global;
                            p.modelnumber=length(zuord);
                            PEAKS{fix}=p;
                            loc_arr(fix)=p.maxlocation;
                            fix=fix+1;
                        end
                        
                    end
                end
                [~,id]=sort(loc_arr);
                PEAKS=PEAKS(id);
                
                obj.peaks(length(PEAKS),1)=peak();
                for i=1:length(PEAKS)
                    if length(PEAKS{i}.peakX_l2r)~=length(PEAKS{i}.fullpeak_l2r)
                        1;
                    end
                    obj.peaks(i)=peak(PEAKS{i}.peakX_l2r,...
                        PEAKS{i}.fullpeak_l2r,...
                        PEAKS{i}.modelnumber,...
                        PEAKS{i}.paramat,...
                        PEAKS{i});
                end
                
                if obj.para.info>=1
                    disp(['### Generate peaks - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'generatePeaks')
            end
        end
        
        %% Post processing peaks - Alt Splitting
        function postprocPeaksSplittingAlt(obj)
            try
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Alt Splitting - Start ###'])
                end
                
                if obj.para.info>=2
                    disp(['=>  Alt. Splitting'])
                end
                                
                % false positive peak position detection
                % a line between two peak position is compared to the signal. if
                % the ration of the maximal orthogonal distance of the signal
                % to the line is larger as ratio limit the smaller peak
                % position is omitted. x values and signal hight is relative.
                ratiolimit=0.1;
                
                %savitzky golay filter settings - order and width (odd)
                sgord=3;
                sgwidth=15;
                
                %signal properties for generation of new peaks
                ma_global=max(obj.profPF.signal);
                
                %%
                act=ones(1,length(obj.peaks));
                while max(act)==1
                    i=find(act==1,1,'first');
                    
                    %get the current peak and its maximum location
                    p=obj.peaks(i);
                    maxloc=p.peakstruct.maxlocation;
                    
                    % fast skip peaks with just one subpeak
                    if p.modelnumber==1
                        act(i)=0;
                        continue;
                    end
                    
                    % fast skip if current or connected peaks are tailing
                    I=obj.getConPeaks(0.01,i);
                    
                    iIsTailing=0;
                    for j=I
                        if obj.peaks(j).istailing
                            iIsTailing=1;
                            break
                        end
                    end
                    if iIsTailing
                        act(i)=0;
                        continue;
                    end
                    
                    % get peak grid/mod
                    modx=p.timegrid;
                    mod=p.signal;
                    
                    % get second derivative by savitzky golay
                    [modsg,modsg1,modsg2]=sgfilt(sgord,sgwidth,mod);
                    
                    % get local minima of sec derivative
                    modsg2scal=-modsg2./max(abs(modsg2));
                    [~,pposidx]=findpeaks(modsg2scal,'MinPeakHeight',0.1);
                    
                    %clear pposidx of false postives
                    cont=1;
                    while cont
                        
                        % check possible peak locations
                        cont=0;
                        for j=1:length(pposidx)-1
                            % get the relevant signal section
                            vec=modsg(pposidx(j):pposidx(j+1));
                            % set minimum to zero
                            vec=vec-min(vec([1 end]));
                            % scale => signal is now in the y-range [0,1]
                            vec=vec./max(vec);
                            
                            %relative position from 0 to 1 by index
                            % signal is now in the x-range [0 1]
                            vecrel=(0:length(vec)-1)./(length(vec)-1);
                            
                            % create a line between start and end point of vec
                            % intersect
                            sup=[0;vec(1)];
                            % direction vector
                            dir=[1;vec(end)-vec(1)];
                            % orthogonal vector
                            ortdir=[-dir(2);dir(1)];
                            % absolute value
                            n3=-(ortdir(2)*sup(2));
                            ndir=norm(ortdir,2);
                            
                            % calc distance between start and end
                            lenpo=norm(dir,2);
                            
                            % calc ratio distance(sig,line) : lenpo
                            distratio=zeros(size(vec));
                            for ii=1:length(distratio)
                                cp=[vecrel(ii);vec(ii)];
                                distratio(ii)=((ortdir'*cp+n3))/ndir/lenpo;
                            end
                            
                            
                            %search for max in ratio (positive values are ommited,
                            %because of wrong curvature). in the case that only
                            %those points are present, this results in ratio=0
                            %and omitting of the smaller peak position
                            
                            ratio=abs(min(distratio));
                            % ratio is always larger than 0
                            
                            if ratio < ratiolimit
                                %ignore the smaller peak position
                                if vec(1)<vec(end)
                                    pposidx(j)=[];
                                    cont=1;
                                    break
                                else
                                    pposidx(j+1)=[];
                                    cont=1;
                                    break
                                end
                            end
                            
                        end
                        
                    end
                    
                    % get splitting locations
                    splitpos=zeros(1,length(pposidx)-1);
                    if length(pposidx)>1
                        for j=1:length(pposidx)-1
                            %split position by 2. deri
                            [~,minidx]=min(modsg2scal(pposidx(j):pposidx(j+1)));
                            %split position by model
                            [~,minidx2]=min(modsg(pposidx(j):pposidx(j+1)));
                            
                            splitpos(j)=pposidx(j)-1+round((minidx+minidx2)/2);
                        end
                    end
                    
                    
                    % if no splitting is needed deactivate peak and continue
                    if isempty(splitpos)
                        act(i)=0;
                        continue;
                    end
                    
                    % split peaks (eventual into multiple peaks)
                    [pArray]=p.splitMulti(modx(splitpos),obj.rawP.timegrid,ma_global);
                    
                    % check cases were splitting results in null peaks
                    % reason: no subpeaks between two splitting positions
                    delidx=[];
                    for j=1:length(pArray)
                        if pArray(j).modelnumber==0
                            delidx(end+1)=j;
                        end
                    end
                    pArray(delidx)=[];
                    
                    if min([pArray(:).modelnumber]) ==0
                        1;
                    end
                    
                    %remove initial peak
                    obj.peaks(i)=[];
                    act(i)=[];
                    
                    %append splitted peaks (deactivated, not to be split again)
                    for j=1:length(pArray)
                        obj.peaks(end+1)=pArray(j);
                        act(end+1)=0;
                    end
                    
                    
                    if obj.para.info>=2
                        disp(['   Splitted peak with max at x=' num2str(maxloc) ' into ' num2str(length(pArray)) ' peaks'])
                    end
                    
                end
                
                
                % resort peaks by maximum location
                if obj.para.info>=2
                    disp('=> Re-Sort')
                end
                
                loc_arr=zeros(1,length(obj.peaks));
                for i=1:length(obj.peaks)
                    loc_arr(i)=obj.peaks(i).peakstruct.maxlocation;
                end
                [~,id]=sort(loc_arr);
                obj.peaks=obj.peaks(id);
                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Alt Splitting - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'postprocPeaksSplittingAlt')
            end
        end
        
        
        %% Post processing peaks - Splitting
        function postprocPeaksSplitting(obj)
            try
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Splitting - Start ###'])
                end
                
                % obj.peaks sind nach center Werten sortiert
                relPeakHi=obj.para.relPeakHi;
                relValley=obj.para.relValley;
                relValleyBoth=obj.para.relValleyBoth;
                
                % para fuer splitting
                sg_wi=40;
                sg_or=3;
                
                
                if obj.para.info>=2
                    disp(['=>  Splitting'])
                end
                
                % Trenne Multi Peaks
                % (Peaks die faelschlicherweise zwei oder mehr offensichtliche
                % Einzelpeaks enthalten. Die Ueberlappung zwischen den Subpeaks
                % war zu gross und fuehrte dazu, dass sie einem Peak zugeordnet
                % wurden.)
                
                ma_global=max(obj.profPF.signal);
                
                
                act=ones(1,length(obj.peaks));
                while max(act)==1
                    
                    for i=1:length(obj.peaks)
                        if (act(i))
                            p=obj.peaks(i);
                            %Lasse peaks aus nur einem Basismodell aus
                            if p.modelnumber==1
                                act(i)=0;
                                continue;
                            end
                            modx=p.timegrid;
                            mod=p.signal;
                            modsg=smooth(mod,sg_wi,'sgolay',sg_or);
                            
                            
                            [ma,id]=max(modsg);
                            
                            %linke seite
                            if id>1
                                lidx=1:id;
                                lidx_rise=find(diff(modsg(lidx)  )<0,1,'last');
                                [tma,tid]=max(mod(1:lidx_rise));
                                lidx_rise_end=tid;

                                issplit=0;
                                
                                if ~isempty(lidx_rise)
                                    %1. Trigger
                                    % split wenn hoehenunterschied bzgl tal gross
                                    % genug und tal tief genug
                                    if mod(lidx_rise_end)-mod(lidx_rise) >relPeakHi*(ma-mod(lidx_rise)) && ...
                                            relValley*ma > mod(lidx_rise)
                                        splitidx=lidx_rise;
                                        issplit=1;
                                    end
                                    
                                    %2. Trigger
                                    if  relValleyBoth*mod(lidx_rise_end) > mod(lidx_rise) && relValleyBoth*ma > mod(lidx_rise)
                                        splitidx=lidx_rise;
                                        issplit=1;
                                    end
                                    
                                end
                                
                                % do split
                                if issplit
                                    [pl,pr]=p.split(modx(splitidx),obj.rawP.timegrid,ma_global);
                                    obj.peaks(i)=[];
                                    act(i)=[];
                                    obj.peaks(end+1)=pl;
                                    obj.peaks(end+1)=pr;
                                    act(end+1)=1;
                                    act(end+1)=1;
                                    
                                    if obj.para.info>=2
                                        disp(['   Splitted peak with max x=' num2str(modx(id)) ' at ' num2str(modx(lidx_rise))])
                                    end
                                    continue
                                end
                                
                            end
                            %rechte seite
                            if id<length(mod)
                                ridx=id:length(mod);
                                ridx_rise=id-1+find(diff(modsg(ridx))>0,1,'first');
                                [tma,tid]=max(mod(ridx_rise:end));
                                ridx_rise_end=ridx_rise-1+tid;
                                issplit=0;
                                
                                if ~isempty(ridx_rise)
                                    %1. Trigger
                                    % split wenn hoehenunterschied bzgl tal gross
                                    % genug und tal tief genug
                                    if mod(ridx_rise_end)-mod(ridx_rise)>relPeakHi*(ma-mod(ridx_rise)) && ...
                                            relValley*ma > mod(ridx_rise)
                                        splitidx=ridx_rise;
                                        issplit=1;
                                    end
                                    
                                    %2. Trigger
                                    if relValleyBoth*mod(ridx_rise_end) > mod(ridx_rise) && relValleyBoth*ma > mod(ridx_rise)
                                        splitidx=ridx_rise;
                                        issplit=1;
                                    end
                                end
                                
                                if issplit
                                    [pl,pr]=p.split(modx(splitidx),obj.rawP.timegrid,ma_global);
                                    obj.peaks(i)=[];
                                    act(i)=[];
                                    obj.peaks(end+1)=pl;
                                    obj.peaks(end+1)=pr;
                                    act(end+1)=1;
                                    act(end+1)=1;
                                    if obj.para.info>=2
                                        disp(['   Splitted peak with max x=' num2str(modx(id)) ' at ' num2str(modx(ridx_rise))])
                                    end
                                    continue
                                end
                            end
                            act(i)=0;
                        end
                    end
                end
                
                
                % Restore sort
                if obj.para.info>=2
                    disp('=> Re-Sort')
                end
                
                loc_arr=zeros(1,length(obj.peaks));
                for i=1:length(obj.peaks)
                    loc_arr(i)=obj.peaks(i).peakstruct.maxlocation;
                end
                [~,id]=sort(loc_arr);
                obj.peaks=obj.peaks(id);
                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Splitting - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'postprocPeaksSplitting')
            end
        end
        
        %% Post processing peaks - Tailing
        function postprocPeaksTailing(obj)
            try                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Tailing - Start ###'])
                end
                % release peaks from tail
                if obj.para.info>=2
                    disp('=> Resolve Tailing')
                end
                
                ma_global=max(obj.profPF.signal);
                
                act=ones(1,length(obj.peaks));
                while max(act)==1
                    for i=1:length(obj.peaks)
                        if (act(i))
                            p=obj.peaks(i);
                            % Skip peaks without possible tailing 
                            % Skip peaks with just one sub model
                            if ~p.istailing
                                act(i)=0;
                                continue;
                            end
                            
                            mod=p.signal;
                            gx=obj.rawP.timegrid;
                            
                            %Get tail length
                            iadd=0;
                            type=[]; %1=tail 2=minipeak
                            typePidx=[];
                            combinedmod=zeros(size(gx));
                            cont=1;
                            while cont
                                cont=0;
                                % is there a next peak?
                                if i+iadd<length(obj.peaks)                                    
                                    if obj.peaks(i+iadd).peakstruct.rbound>obj.peaks(i+iadd+1).peakstruct.lbound
                                        
                                        mod1=obj.peaks(i+iadd).getmod(gx);
            
                                        if iadd==0
                                            combinedmod=mod1;
                                        end
                                        
                                        mod2=obj.peaks(i+iadd+1).getmod(gx);
                                        I_t=trapz(gx,combinedmod);
                                        I_t2=trapz(gx,mod2);
                                        I_inter=trapz(gx,min(combinedmod,mod2));
                                        
                                        % is next peak overlapping and not
                                        % completely contained and not
                                        % larger than current peak
                                        [mcmax,~]=max(combinedmod);
                                        [m2,i2]=max(mod2);
                                        if I_inter/I_t2>1E-3 && ...         % overlap to small
                                                I_inter/I_t2<=1 && ...      % overlap to large
                                                0.3*mcmax >= m2             % new peak is smaller than current main peak
                                            [m1,i1]=max(mod1);
                                            iadd=iadd+1;
                                            modc=mod1+mod2;
                                            mc1=modc(i1);
                                            mc2=modc(i2);
                                            mic=min(modc(min(i1,i2):max(i1,i2)));
                                            
                                            combinedmod=combinedmod+mod2;
                                            if mc2/mic>1.05
                                                %minipeak
                                                type(end+1)=2;
                                            else
                                                %tail
                                                type(end+1)=1;
                                            end
                                            typePidx(end+1)=i+iadd;
                                            cont=1;
                                        end
                                    end
                                end
                            end
                            
                            % cut mini peaks at end
                            if ~isempty(type)
                                while type(end)==2
                                    type(end)=[];
                                    typePidx(end)=[];
                                    if isempty(type)
                                        break
                                    end
                                end
                            end
                            
                            % split into peak in tail and mini peaks on tail
                            
                            %tailing part included
                            if isempty(type)
                                act(i)=0;
                                continue
                            end
                            % sort 
                            if max(type==1) && max(type==2)
                                % start with main peak
                                % p = current main peak
                                mainp=p;
                                clear minip
                                minip(1)=peak();
                                triggered=0;
                                
                                for ii=1:length(type)
                                    if ii>1
                                        if type(ii)==2 && type(ii-1)==2
                                            % break for two subsequent peaks
                                            type(ii:end)=[];
                                            typePidx(ii:end)=[];
                                            break
                                        end
                                    end
                                    if type(ii)==1
                                        mainp(end+1)=obj.peaks(i+ii);
                                    else
                                        if length(minip)==1 && ~triggered
                                            triggered=1;
                                            minip(1)=obj.peaks(i+ii);
                                        else
                                            minip(end+1)=obj.peaks(i+ii);
                                        end
                                    end
                                end
                                % Determine peak with tail without mini
                                % peaks
                                % -> modflat
                                modflat=zeros(size(gx));
                                paraflat=[];
                                for ii=1:length(mainp)
                                    modflat=modflat+mainp(ii).getmod(gx);
                                    paraflat=[paraflat mainp(ii).parameters];
                                end
                                
                                for ii=1:length(minip)
                                    addidx=0;
                                    redu=0.10;
                                    lb=minip(ii).peakstruct.lbound;
                                    rb=minip(ii).peakstruct.rbound;
                                    leng=rb-lb;
                                    lb=lb+leng*redu;
                                    rb=rb-leng*redu;
                                    idxlb=find(gx>lb,1,'first');
                                    idxrb=find(gx>=rb,1,'first');
                                    idxlb=max(idxlb-addidx,1);
                                    idxrb=min(idxrb+addidx,length(gx));
                                    ly=modflat(idxlb);
                                    ry=modflat(idxrb);
                                    gapy=linspace(ly,ry,idxrb-idxlb+1);
                                    
                                    gapmod=zeros(size(gx));
                                    gapmod(idxlb:idxrb)=gapy;
                                    modflat=max(modflat,gapmod);
                                end
                                
                                %Build result peaks
                                %main peak
                                param=paraflat;
                                
                                modeln=size(param,2);
                                
                                
                                % determine peak with correction
                                % then update properties
                                fullpeak=zeros(length(gx),1);
                                for jj=1:modeln
                                    fullpeak=fullpeak+mod_ga(gx,param(:,jj));
                                end
                                
                                [ma,id]=max(fullpeak);
                                left=find(fullpeak>ma/300,1,'first');
                                right=find(fullpeak>ma/300,1,'last');
                                if isempty(left)
                                    left=1;
                                end
                                if isempty(right)
                                    right=length(fullpeak);
                                end
                                inte=trapz(gx,fullpeak);
                                if obj.para.info>=2
                                    disp( '=> Tailing main peak        ');
                                    disp(['     Left bound:        ' num2str(gx(left))]);
                                    disp(['     Peak max location: ' num2str(gx(id))]);
                                    disp(['     Right bound:       ' num2str(gx(right))]);
                                    disp(['     Peak height:       ' num2str(ma)]);
                                    disp(['     Peak rel. height:  ' num2str(ma/ma_global)]);
                                    disp(['     Model number:      ' num2str(modeln)]);
                                end
                                pp.peakX_l2r=gx(left:right);
                                pp.fullpeak_l2r=fullpeak(left:right);
                                pp.paramat=param;
                                pp.lbound=gx(left);
                                pp.rbound=gx(right);
                                pp.maxlocation=gx(id);
                                pp.height=ma;
                                pp.relheight=ma/ma_global;
                                pp.modelnumber=modeln;
                                
                                mainpeak=peak(pp.peakX_l2r,pp.fullpeak_l2r,modeln,param,pp);
                                mainpeak.correction=interp1(gx,modflat-mainpeak.getmod(gx),pp.peakX_l2r,'linear',0);
                                mainpeak.updatePropAfterCorr(gx,ma_global);
                                
                                %correction of mini peaks
                                for jj=1:length(minip)
                                    cgr=minip(jj).peakstruct.peakX_l2r;
                                    cgrF=pp.peakX_l2r;
                                    ccorF=mainpeak.correction;
                                    localcor=interp1(cgrF,ccorF,cgr,'linear',0);
                                    
                                    
                                    minip(jj).correction=-localcor;
                                    minip(jj).updatePropAfterCorr(gx,ma_global);
                                end
                                
                                %update global peaks
                                del=[i typePidx(:)'];
                                act(del)=[];
                                obj.peaks(del)=[];
                                obj.peaks(end+1)=mainpeak;
                                act(end+1)=0;
                                
                                
                                for jj=1:length(minip)
                                    obj.peaks(end+1)=minip(jj);
                                    act(end+1)=0;
                                end
                            else
                                act(i)=0;
                            end
                            break
                        end
                    end
                end
                
                % resort peaks according to their location
                if obj.para.info>=2
                    disp('=> Re-Sort')
                end
                
                loc_arr=zeros(1,length(obj.peaks));
                for i=1:length(obj.peaks)
                    loc_arr(i)=obj.peaks(i).peakstruct.maxlocation;
                end
                [~,id]=sort(loc_arr);
                obj.peaks=obj.peaks(id);
                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Tailing - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'postprocPeaksTailing')
            end
        end
        
        
        %% Post processing peaks - Adjust Tails
        function postprocPeaksAdjustTails(obj)
            try
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Adjust Tailing - Start ###'])
                end
                
                ma_global=max(obj.profPF.signal);
                fullgx=obj.profPF.timegrid;
                
                % find possible tailing peaks
                Tcluster={};
                for i=1:length(obj.peaks)
                    if obj.peaks(i).istailingDetail(obj.profPF,obj.rescl.noiselvl)
                        Tcluster{end+1}=i;
                    end
                end
                if obj.para.info>=1
                    if ~isempty(Tcluster)
                        disp('=> found tailing peaks:')
                        disp(Tcluster)
                    end
                end
                
                %find overlapping predecessor peaks
                for i=1:length(Tcluster)
                    %check overlap with predecessor
                    cidx=Tcluster{i};
                    while 1
                        % no previous peaks
                        if cidx==1
                            break
                        end
                        
                        mx1=obj.peaks(cidx-1).peakstruct.maxlocation;
                        mx2=obj.peaks(cidx).peakstruct.maxlocation;
                        rb1=obj.peaks(cidx-1).peakstruct.rbound;
                        lb2=obj.peaks(cidx).peakstruct.lbound;
                        midx1=find(fullgx>=mx1,1,'first');
                        midx2=find(fullgx<=mx2,1,'last');
                        
                        %ensure that there's at least one intermediate point
                        if midx1+1<midx2
                            %get tops of both peaks
                            mh1=obj.profPF.signal(midx1);
                            mh2=obj.profPF.signal(midx2);
                            
                            %get minimum between tops
                            interSig=obj.profPF.signal(midx1:midx2);
                            minHi=min(interSig);
                            
                            % check if valley reaches 0 and predecessor should
                            % not be to large
                            if minHi>=0.05*min(mh1,mh2) && mh1<=2*mh2  && rb1>lb2
                                isoverlap=1;
                            else
                                isoverlap=0;
                            end
                        else
                            %no intermediate point
                            isoverlap=0;
                        end
                        
                        if isoverlap
                            cidx=cidx-1;
                        else
                            break
                        end
                    end
                    
                    %overlapping peaks exist
                    %add to Tcluster
                    if cidx<Tcluster{i}
                        Tcluster{i}=cidx:Tcluster{i};
                    end
                end
                
                % remove single tailing peaks
                delidx=[];
                for i=1:length(Tcluster)
                    if length(Tcluster{i})<=1
                        delidx(end+1)=i;
                    end
                end
                Tcluster(delidx)=[];
                
                % print
                if obj.para.info>=1
                    if ~isempty(Tcluster)
                        disp('=> overlapping predecessors:')
                        for i=1:length(Tcluster)
                            disp(['<> Cluster ' num2str(i) ':'])
                            disp(Tcluster{i})
                        end
                    end
                end
                
                % modify cluster peaks
                alreadyAdapted=zeros(1,length(obj.peaks));
                
                for i=1:length(Tcluster)
                    %get last peak in a cluster
                    % leadP is a reference
                    % changes are applied to obj.peaks(.) as well
                    leadP=obj.peaks(Tcluster{i}(end));
                    sig=leadP.getmod(leadP.timegrid);
                    
                    %extract tailing peak shape information
                    scaleTG  = linspace(0,1,length(leadP.timegrid));
                    scaleSig = sig./max(sig);
                    
                    % get point with where a tangent with slope tSlope touches signal
                    % on the right side of the maximum
                    
                    [~,maIdx]=max(scaleSig);
                    
                    %get slope of line from top to end of tail
                    tSlope=diff(scaleSig([maIdx end]))/diff(scaleTG([maIdx end]));
                    
                    %calc y-intersec of lines through(scaleTG(.),scaleSig(.))
                    %with slope tSlope
                    % (todo: calc only relevant points)
                    inter=scaleSig(:)-tSlope*scaleTG(:);
                    inter(1:maIdx-1)=inf;
                    [~,tIdx]=min(inter);
                    
                    tPoint=[scaleTG(tIdx),scaleSig(tIdx)];
                    
                    %ratio of tail length to signal height
                    tailratio=(leadP.timegrid(end)-leadP.timegrid(tIdx))/max(sig);
                    
                    % apply shape to previous peaks in cluster
                    %for j=length(Tcluster{i})-1:-1:1
                    for j=length(Tcluster{i})-1
                        if j<=length(Tcluster{i})-2
                            break
                        end
                        
                        % cPeak is a reference
                        % changes are applied to obj.peaks(.) as well
                        cPeak=obj.peaks(Tcluster{i}(j));
                        cSig=cPeak.getmod(cPeak.timegrid);
                        [maCPeak,maIdxCPeak]=max(cSig);
                        tStartIdx=find(cPeak.timegrid >= cPeak.timegrid(maIdxCPeak) & cSig <= maCPeak*tPoint(2),1,'first');
                        if isempty(tStartIdx)
                            break
                        end
                        
                        % generate tail
                        tailmodx=0:0.02:1;
                        tailmods=interp1([0 0.1 0.2 1],[1 0.5 0.3 0],tailmodx,'pchip');
                        
                        %start in cpeak max location
                        tstart=cPeak.timegrid(maIdxCPeak);
                        %end is correlated to peak height
                        tlengthInf=tailratio*maCPeak;
                        %maximum length is limited by leadP rbound
                        tendInf=tstart+tlengthInf;
                        tendInf=min(tendInf,leadP.peakstruct.rbound);
                        tlengthRes=tendInf-tstart;
                        
                        % safty margin (substraction of correction for leadP will
                        % push rbound to the left)
                        tlength=0.96*tlengthRes;
                        
                        %scaling/shifting grid
                        tailmodx=tailmodx*tlength;
                        tailmodx=tailmodx+tstart;
                        
                        %scaling tail signal
                        tail=tailmods*0.95*cSig(tStartIdx);
                        
                        %get tail on full grid
                        cTailFull=interp1(tailmodx,tail,fullgx,'linear',0);
                        
                        %get cPeak signal on full grid
                        cSigOri=cPeak.getmod(fullgx);
                        
                        %get cPeak signal and tail on full grid
                        cSigFull=max(cTailFull,cSigOri);
                        
                        %ensure correction is below cPeak and leadP
                        cSigFull=min(cSigFull,cSigOri+leadP.getmod(fullgx));
                        
                        %get correction on full grid
                        cSigDiff=cSigFull-cSigOri;
                        
                        % insert corrections
                        cPeak.addCorrection(fullgx,ma_global,cSigDiff);
                        leadP.addCorrection(fullgx,ma_global,-cSigDiff);
                    end
                    
                end
                
                obj.updateModelAndError;
                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Adjust Tailing - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'postprocPeaksAdjustTails')
            end
        end
        
        
        %% Post process local fitting
        function obj=postLocalFit(obj,startidx)
            try
                I=obj.getConPeaks(0.01,startidx);
                if obj.para.info>=2
                    disp(['=> Conncected Peaks ' num2str(length(I))])
                end
                gx=obj.profPF.timegrid;
                sig=obj.profPF.signal;
                
                % get bounds
                lb=obj.peaks(startidx).peakstruct.lbound;
                rb=obj.peaks(startidx).peakstruct.rbound;
                for i=I
                    %subs correction
                    sig=sig-interp1(obj.peaks(i).timegrid,obj.peaks(i).correction,gx,'linear',0);
                    %bounds
                    lbn=obj.peaks(i).peakstruct.lbound;
                    rbn=obj.peaks(i).peakstruct.rbound;
                    if lbn<=lb
                        lb=lbn;
                    end
                    if rbn>=rb
                        rb=rbn;
                    end
                end
                idxL=find(gx>=lb,1,'first');
                idxR=find(gx<=rb,1,'last');
                
                %get parameters
                param=[];
                for i=I
                    paraadd=obj.peaks(i).parameters;
                    
                    %parameter merge
                    param=[param paraadd];
                    
                end

                %Fit START
                currtemp=sig(idxL:idxR);
                scaltemp=max(abs(currtemp));
                datf=(1/scaltemp)*currtemp;
                
                
                %bounds
                param(4,:)=(1/scaltemp)*param(4,:);
                paramLB=param;
                paramUB=param;
                %center
                paramLB(1,:)=param(1,:)-4*param(2,:); %center- 4 sigmaL
                paramUB(1,:)=param(1,:)+4*param(3,:); %center- 4 sigmaR
                %sigma left
                paramLB(2,:)=0.2*paramLB(2,:);
                paramUB(2,:)=3*param(2,:);
                %sigma right
                paramLB(3,:)=0.2*paramLB(3,:);
                paramUB(3,:)=3*param(3,:);
                %height
                %             paramLB(4,:)=0.1*param(4,:);
                %             paramUB(4,:)=1.5*param(4,:);
                paramLB(4,:)=0.6*param(4,:);
                paramUB(4,:)=1.5*param(4,:);
                %linearity left/right
                paramLB(5:6,:)=0;
                paramUB(5:6,:)=1;
                
                
                %Change of initial param
                param(4,:)=1*param(4,:);
                
                gxtemp=gx(idxL:idxR);
                
                %exclude ranges
                exIdxRange=obj.getExcludeIndices(gxtemp);
                exIdxRange=exIdxRange(:)';
                actx=true(1,length(gxtemp)) & exIdxRange;
                actx_w=true(1,length(gxtemp)) & exIdxRange;
                
                
                
                
                we=[1 20];
                fun2=@(p) tar_opt(p,gxtemp,datf,size(param,2),6,we,[],[],actx,actx_w,0);
                options=optimoptions(@lsqnonlin,...
                    'Display','off',...
                    'MaxIter',50,...
                    'MaxFunEvals',1000,...
                    'TolFun',1e-10,...
                    'TolX',1e-10);
                pOpt=lsqnonlin(fun2,param,paramLB,paramUB,options);
                paramOpt=reshape(pOpt,size(param,1),size(param,2));
                paramOpt(4,:)=scaltemp*paramOpt(4,:);
                
                id=0;
                ma_global=max(obj.profPF.signal);
                ii=1;
                for i=I
                    obj.peaks(i).parameters=paramOpt(:,id+1:id+size(obj.peaks(i).parameters,2));
                    obj.peaks(i).updatePropAfterCorr(gx,ma_global);
                    id=id+size(obj.peaks(i).parameters,2);
                    ii=ii+1;
                end
                
            catch ex
                if ex.identifier == 'optimlib:lsqncommon:ProblemNotHandled'
                    disp('<> fitting problem in postLocalFit: max error skipped')
%                     obj.excludeRange=[obj.excludeRange; idxL idxR];
                else
                	obj.errorLog(ex,'postLocalFit')
                end
            end
        end
        
        %% filter small and narrow peaks
        function filterPeaks(obj,nlvlfac,ind)
            try
                if obj.para.info>=1
                    disp(['### Filter peaks - Start ###'])
                end
                oldanz=length(obj.peaks);
                del=zeros(length(obj.peaks),1);
                idxv=1:length(obj.peaks);
                idxv=idxv(:);
                for i=1:length(obj.peaks)
                    %                 if obj.peaks(i).peakstruct.maxlocation>102
                    %                    1;
                    %                 end
                    
                    % filter narrow peaks
                    if length(obj.peaks(i).peakstruct.peakX_l2r)<ind
                        del(i)=1;
                        if obj.para.info>=2
                            disp(['=> filter narrow peak at: ' num2str(obj.peaks(i).peakstruct.maxlocation)])
                        end
                    end
                    % filter small peaks
                    %                 lb_idx=find(obj.profPF.timegrid>=obj.peaks(i).peakstruct.lbound,1,'first');
                    %                 lb_sig=obj.profPF.signal(lb_idx);
                    %                 rb_idx=find(obj.profPF.timegrid<=obj.peaks(i).peakstruct.rbound,1,'last');
                    %                 rb_sig=obj.profPF.signal(rb_idx);
                    %                 m_sig=(lb_sig+rb_sig)/2;
                    %                 if obj.peaks(i).peakstruct.height-m_sig<=nlvlfac*obj.rescl.noiselvl
                    if obj.peaks(i).peakstruct.height<=nlvlfac*obj.rescl.noiselvl
                        del(i)=1;
                        if obj.para.info>=2
                            disp(['=> filtered small peak at: ' num2str(obj.peaks(i).peakstruct.maxlocation)])
                        end
                    end
                end
                
                del=logical(del);
                modsubs=zeros(size(obj.profPF.timegrid));
                for i=1:length(obj.peaks)
                    if del(i)
                        modsubs=modsubs+obj.peaks(i).getmod(obj.profPF.timegrid);
                    end
                end
                obj.modePF.signal=obj.modePF.signal-modsubs;
                
                obj.peaks(logical(del))=[];
                length(obj.peaks)
                if obj.para.info>=1
                    
                    newanz=length(obj.peaks);
                    disp(['=> Filtered ' num2str(oldanz-newanz) ' of ' num2str(oldanz) ' peaks'])
                    disp(['=> New number of peaks: ' num2str(newanz)])
                    disp(['### Filter peaks - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'filterPeaks')
            end
        end
        
        
        %% Get peak indices in range [lloc , rloc]
        function I=getPeaks(obj,lloc,rloc)
            I=[];
            if ~isempty(obj.peaks)
                for i=1:length(obj.peaks)
                    pp=obj.peaks(i).peakstruct;
                    if lloc<=pp.maxlocation && pp.maxlocation<=rloc
                        I(end+1)=i;
                    end
                end
            else
                disp('No peaks in range.')
            end
        end
        
        %% Get closest peak to gx(idx)
        function minIdx=getClosestPeak2IDX(obj,idx)
            D=inf(length(obj.peaks),1);
            gx=obj.profPF.timegrid;
            gxidx=gx(idx);
            minIdx=[];
            if ~isempty(obj.peaks)
                for i=1:length(obj.peaks)
                    pp=obj.peaks(i).peakstruct;
                    D(i)=abs(pp.maxlocation-gxidx);
                end
                [~,minIdx]=min(D);
            else
                disp('No peaks calculated.')
            end
        end
        
        %% Get closest peak to gx(idx)
        function minIdx=getClosestPeak2GX(obj,gxidx)
            D=inf(length(obj.peaks),1);
            if ~isempty(obj.peaks)
                for i=1:length(obj.peaks)
                    pp=obj.peaks(i).peakstruct;
                    D(i)=abs(pp.maxlocation-gxidx);
                end
                [~,minIdx]=min(D);
            else
                disp('No peaks calculated.')
            end
        end
        
        %% Get connected p. with >= (100*minover) % overlap starting at p. startidx
        function I=getConPeaks(obj,minover,startidx)
            gx=obj.modePF.timegrid;
            I=startidx;
            countlimit=8;
            
            combinedmod=obj.peaks(startidx).getmod(gx);
            %right side
            
            count=0;
            for ir=startidx+1:length(obj.peaks)
                mod2=obj.peaks(ir).getmod(gx);
                I_t2=trapz(gx,mod2);
                I_inter=trapz(gx,min(combinedmod,mod2));
                if I_inter/I_t2>=minover
                    I(end+1)=ir;
                    combinedmod=combinedmod+mod2;
                else
                    count=count+1;
                    if count>countlimit
                        break
                    end
                end
            end
            
            
            %left side
            count=0;
            for il=startidx-1:-1:1
                mod2=obj.peaks(il).getmod(gx);
                I_t2=trapz(gx,mod2);
                I_inter=trapz(gx,min(combinedmod,mod2));
                if I_inter/I_t2>=minover
                    I(end+1)=il;
                    combinedmod=combinedmod+mod2;
                else
                    count=count+1;
                    if count>countlimit
                        break
                    end
                end
            end
            
            
            I=sort(unique(I));
        end
        
        %% Write peaks
        function writePeaks(obj,fpath)
            if nargin==1
                fpath=[obj.para.tempf filesep 'peak_res.csv'];
            end
            fid = fopen(fpath,'w');
            fprintf(fid,'%12s, ',...
                'index',...
                'max location',...
                'left bound',...
                'right bound',...
                'height',...
                'rel. height',...
                'area',...
                'rel area',...
                'sig area',...
                'rel sig area',...
                'agi area',...
                'rel agi area');
            fprintf(fid,'%12s','subpeak number');
            fprintf(fid,'\n');
            for i=1:length(obj.peaks)
                
                fprintf(fid,'%12f, ',i ,obj.peaks(i).peakstruct.maxlocation, ...
                    obj.peaks(i).peakstruct.lbound, ...
                    obj.peaks(i).peakstruct.rbound, ...
                    obj.peaks(i).peakstruct.height, ...
                    obj.peaks(i).peakstruct.relheight,...
                    obj.peaks(i).intinfo.abs, ...
                    obj.peaks(i).intinfo.rel, ...
                    obj.peaks(i).intinfo.sigAbs, ...
                    obj.peaks(i).intinfo.sigRel, ...
                    obj.peaks(i).intinfo.agiAbs, ...
                    obj.peaks(i).intinfo.agiRel);
                fprintf(fid,'%12f',obj.peaks(i).peakstruct.modelnumber);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        %% evaluate peaks for overlap and fit error
        function evalPeaks(obj,idx)
            try
                if nargin==1
                    range=1:length(obj.peaks);
                else
                    range=idx;
                end
                
                for i=range
                    pp=obj.peaks(i);
                    cmod=pp.getmod(obj.profPF.timegrid);
                    cmod=cmod(:);
                    
                    intC=trapz(obj.profPF.timegrid, cmod);
                    plo=0;
                    pro=0;
                    %left peak overlap
                    if i>1
                        lmod=obj.peaks(i-1).getmod(obj.profPF.timegrid);
                        intL=trapz(obj.profPF.timegrid, min(lmod,cmod));
                        plo=intL/(intC);
                    end
                    %right peak overlap
                    if i<length(obj.peaks)
                        rmod=obj.peaks(i+1).getmod(obj.profPF.timegrid);
                        intR=trapz(obj.profPF.timegrid, min(rmod,cmod));
                        pro=intR/(intC);
                    end
                    obj.peaks(i).overlap=max(plo,pro);
                    
                    ce=pp.peakstruct.maxlocation;
                    wi=(pp.peakstruct.rbound-pp.peakstruct.lbound)/4;
                    leval=find(obj.profPF.timegrid>=ce-wi,1,'first');
                    reval=find(obj.profPF.timegrid<=ce+wi,1,'last');
                    evalsig=obj.profPF.signal(leval:reval);
                    evalmod=obj.modePF.signal(leval:reval);
                    evalres=evalsig-evalmod;
                    
                    obj.peaks(i).fit=sqrt((evalres'*evalres)/(evalsig'*evalsig));
                end
                
            catch ex
                obj.errorLog(ex,'evalPeaks')
            end
        end
        
        
        %% start PostProc Analyzer
        function startAnalyzer(obj)
            RESanalysis2(obj);
        end
        
        %% Plotting
        function plotPeaks(obj)
            figure
            plot(obj.profPF,'k',2);
            hold on
            plot(obj.modePF,'b',1);
            axis tight
            yy=ylim;
            ylim([min(yy(2),0) yy(2)])
            
            cc =[0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
            
            id=1;
            for i=1:length(obj.peaks)
                plot(obj.peaks(i),cc(id,:),2);
                id=id+1;
                if id==8
                    id=1;
                end
                
            end
        end
        
        function plotBL(obj)
            figure
            plot(obj.rawP,'b');
            hold on
            plot(obj.baseP,'k');
            axis tight
            yy=ylim;
            ylim([min(yy(2),0) yy(2)])
            legend('raw data','base. cor. data')
        end
        
        
        function plot(obj,cfig)
            
            marktailing=1;
            if isnan(obj.progress(2))
                disp('nothing to plot. calculate something')
                return
            end
            
            if nargin==2
                if ~ishandle(cfig)
                    error('second argument has to be a handle')
                end
            end
            
            if nargin==1
                cfig=figure;
            end
            
            
            if nargin==1
                set(gcf,'position',[100 100 1400 800])
            end
            if ~isnan(obj.progress(2))
                if ~isnan(obj.progress(3))
                    ax1=subplot(3,2,1);
                end
                
                plot(obj.rawP,'b');
                hold on
                plot(obj.baseP,'k');
                axis tight
                yy=ylim;
                ylim([min(yy(2),0) yy(2)])
                legend('raw data','base. cor. data')
            end
            
            
            
            if ~isnan(obj.progress(3))
                ax2=subplot(3,2,2);
                plot(obj.baseP,'k');
                hold on
                axis tight
                yy=ylim;
                ylim([min(yy(2),0) yy(2)])
                
                for i=1:length(obj.clusterP)
                    plot(obj.basePF,'b');
                    lx=obj.clusterP(i).timegrid(1);
                    rx=obj.clusterP(i).timegrid(end);
                    fill([lx rx rx lx],yy([1 1 2 2]),'b','facecolor','b','facealpha',0.1);
                end
                legend('base. cor. data','local baseline','cluster')
            end
            
            if ~isnan(obj.progress(4))
                ax3=subplot(3,3,4:9);
                plot(obj.profPF,'k',2);
                hold on
                plot(obj.modePF,'b',1);
                axis tight
                yy=ylim;
                ylim([min(yy(2),0) yy(2)])
                
                
                
                title(['SSQ-Error: ' num2str(obj.modelErrorNoN) ' (' num2str(obj.modelError) ...
                    '); Explained Signal Area: ' num2str(num2str(obj.coveredPeakAreaNoN)) ' (' num2str(num2str(obj.coveredPeakArea))...
                    '); NL: ' num2str(obj.rescl.noiselvl)]);
                
                cc =[0    0.4470    0.7410
                    0.8500    0.3250    0.0980
                    0.9290    0.6940    0.1250
                    0.4940    0.1840    0.5560
                    0.4660    0.6740    0.1880
                    0.3010    0.7450    0.9330
                    0.6350    0.0780    0.1840];
                
                id=1;
                for i=1:length(obj.peaks)
                    plot(obj.peaks(i),cc(id,:),2);
                    id=id+1;
                    if id==8
                        id=1;
                    end
                    if marktailing
                        if obj.peaks(i).istailingDetail(obj.profPF, obj.rescl.noiselvl)
                            ce=obj.peaks(i).peakstruct.maxlocation;
                            hi=obj.peaks(i).peakstruct.height;
                            plot(ce,hi*1.05,'or','linewidth',2);
                        end
                    end
                end
                legend('base. cor. data (loca+glob)','model fit')
            end
            
            if exist('ax3')
                linkaxes([ax1 ax2 ax3])
            else
                if exist('ax2')
                    linkaxes([ax1 ax2])
                end
            end
        end
        
        %% Plotting (just main subplot)
        function f=plotSingle(obj,show)
            marktailing=1;
            if isnan(obj.progress(2))
                disp('nothing to plot. calculate something')
                return
            end
            if nargin==1
                show=1;
            end
            if show
                f=figure;
            else
                f=figure('visible','off');
            end
            set(gcf,'position',[100 100 1400 800])
            
            
            if ~isnan(obj.progress(4))
                plot(obj.profPF,'k',2)
                hold on
                plot(obj.modePF,'b',1)
                axis tight
                yy=ylim;
                ylim([min(yy(2),0) yy(2)])
                
                
                
                title(['SSQ-Error: ' num2str(obj.modelErrorNoN) ' (' num2str(obj.modelError) ...
                    '); Explained Signal Area: ' num2str(num2str(obj.coveredPeakAreaNoN)) ' (' num2str(num2str(obj.coveredPeakArea))...
                    ')']);
                
                cc =[0    0.4470    0.7410
                    0.8500    0.3250    0.0980
                    0.9290    0.6940    0.1250
                    0.4940    0.1840    0.5560
                    0.4660    0.6740    0.1880
                    0.3010    0.7450    0.9330
                    0.6350    0.0780    0.1840];
                
                id=1;
                for i=1:length(obj.peaks)
                    plot(obj.peaks(i),cc(id,:),2);
                    id=id+1;
                    if id==8
                        id=1;
                    end
                    if marktailing
                        if obj.peaks(i).istailing
                            ce=obj.peaks(i).peakstruct.maxlocation;
                            hi=obj.peaks(i).peakstruct.height;
                            plot(ce,hi*1.05,'or','linewidth',2)
                        end
                    end
                end
                legend('base. cor. data (loca+glob)','model fit')
            end
        end
        
        
        %% Plot Reference Comparison
        function [ax1,ax2]=plotRefComp(obj,isAll)
            if nargin==1
                isAll=0;
            end
            %             if isAll
            %                 isArea=0;
            %                 disp('Warning: Area comparison is turned off, if all Peaks are shown')
            %             end
            if ~isnan(obj.progress(4)) && obj.para.isREF
                figure
                ax1=subplot(3,1,1:2);
                plot(obj.profPF,'k');
                hold on
                plot(obj.modePF,'b',2);
                axis tight
                yy=ylim;
                ylim([min(yy(2),0) yy(2)])
                
                
                
                title(['SSQ-Error: ' num2str(obj.modelErrorNoN) ' (' num2str(obj.modelError) ...
                    '); Explained Signal Area: ' num2str(num2str(obj.coveredPeakAreaNoN)) ' (' num2str(num2str(obj.coveredPeakArea))...
                    ')']);
                
                cc =[0    0.4470    0.7410
                    0.8500    0.3250    0.0980
                    0.9290    0.6940    0.1250
                    0.4940    0.1840    0.5560
                    0.4660    0.6740    0.1880
                    0.3010    0.7450    0.9330
                    0.6350    0.0780    0.1840];
                
                
                
                dataR=obj.para.dataR;
                if size(dataR,2)==7
                    isnew=0;
                elseif size(dataR,2)==2
                    isnew=1;
                else
                    error('reference format unknown')
                end
                
                if ~isnew
                    ref_nr      = size(dataR,1);
                    ref_pos     = dataR(:,2);
                    ref_hi      = dataR(:,4);
                    ref_int     = dataR(:,3);
                    ref_intrel     = dataR(:,6);
                else
                    ref_nr      = size(dataR,1);
                    ref_pos     = dataR(:,1);
                    ref_hi      = ones(ref_nr,1); % not available
                    ref_int     = dataR(:,2);
                    ref_intrel  = dataR(:,2)./sum(dataR(:,2));
                end
                [~,idx]=sort(ref_pos);
                ref_pos     = ref_pos(idx);
                ref_hi      = ref_hi(idx);
                ref_int     = ref_int(idx);
                ref_intrel  = ref_intrel(idx);
                
                
                for i=1:ref_nr
                    plot(ref_pos(i)*[1 1],[0 ref_hi(i)],'bo-','linewidth',2,'markersize',13);
                end
                
                %plot extract
                for j=1:length(obj.peaks)
                    if isAll
                        param=obj.peaks(j).parameters;
                        for i=1:size(param,2)
                            p=param(:,i);
                            p_hi=p(4);
                            p_pos=p(1);
                            plot(p_pos*[1 1],[0 p_hi],'rx-');
                        end
                    else
                        ps=obj.peaks(j).peakstruct;
                        hi=ps.height;
                        plot(ps.maxlocation*[1 1],[0 hi],'rx-','linewidth',2,'markersize',13);
                    end
                end
                
                %Vorbereitung 2.Plot
                %Binning (Mitten zwischen Peakmax Ref-Data)
                binidx=1;
                gx=obj.profPF.timegrid;
                for i=1:ref_nr-1
                    binx=mean(ref_pos(i:i+1));
                    binidx(end+1)=find(gx>=binx,1,'first');
                end
                binidx(end+1)=length(gx);
                
                binres=zeros(size(binidx));
                
                for i=1:length(obj.peaks)
                    cp=obj.peaks(i);
                    ml=cp.peakstruct.maxlocation;
                    if ml>60.4 && ml<60.55
                        disp('a')
                    end
                    id=find(gx(binidx)>=ml,1,'first')-1;
                    binres(id)=binres(id)+1;
                end
                
                ax2=subplot(3,1,3);
                plot(gx(binidx),binres-1,'b','Linewidth',2);
                hold on
                plot(gx([1 end]),[0,0],'r--','Linewidth',2);
                yy=ylim;
                ylim([-1 yy(2)])
                
                linkaxes([ax1 ax2],'x');
                axes(ax1)
                
                %legend('base. cor. data (loca+glob)','model fit')
            else
                warning('Complete computation first OR no reference loaded')
            end
            
        end
        
        %% Illustrate noise level
        function printNoiseLvl(obj)
            noiselvl=obj.rescl.noiselvl;
            
            figure;
            set(gcf,'position',[100 100 1400 800])
            plot(obj.profPF,'k',1)
            hold on
            %             plot(obj.modePF,'b',1)
            xx=xlim();
            plot(xx,noiselvl*[1 1],'r')
            axis tight
            ylim([0 15*noiselvl]);
        end
        
        %% Get Error (exclusive exRanges)
        function E=getError(obj,tol)
            if nargin==1
                tol=0.05;
            end
            sig=obj.profPF.signal;
            mod=obj.modePF.signal;
            E=abs(mod-sig)./max(abs(mod),tol*max(abs(mod)));
            
            %remove exclude ranges (no arg=full range)
            exIdxRange=obj.getExcludeIndices;
            E(~exIdxRange)=0;
        end
        
        %% Illustrate Error (exclusive exRanges)
        function plotError(obj)
            x=obj.profPF.timegrid;
            
            figure;
            set(gcf,'position',[100 100 1400 800])
            
            ax1=subplot(2,1,1);
            plot(obj.profPF,'k',1)
            hold on
            plot(obj.modePF,'b',1)
            axis tight
            
            ax2=subplot(2,1,2);
            plot(x,obj.getError)
            linkaxes([ax1 ax2],'x');
        end
        
        %% Print Peak info (optional: in a certain range)
        function dispPeaks(obj,lloc,rloc)
            % show peak info
            % call ga.dispPeak or ga.dispPeak(lloc,rloc)
            % to show all peaks or those in the interval [lloc,rloc]
            
            Irestrict=1:length(obj.peaks);
            
            if nargin==1
                lloc=-inf;
                rloc=inf;
            elseif nargin==2
                Irestrict=lloc;
                lloc=-inf;
                rloc=inf;
            end
            
            
            I=obj.getPeaks(lloc,rloc);
            
            if ~isempty(I)
                for i=I
                    if max(Irestrict==i)
                        pp=obj.peaks(i).peakstruct;
                        
                        disp(['peak #' num2str(i)]);
                        disp(['    Left bound:        ' num2str(pp.lbound)]);
                        disp(['    Peak max location: ' num2str(pp.maxlocation)]);
                        disp(['    Right bound:       ' num2str(pp.rbound)]);
                        disp(['    Peak height:       ' num2str(pp.height)]);
                        disp(['    Peak rel. height:  ' num2str(pp.relheight)]);
                        disp(['    Model number:      ' num2str(pp.modelnumber)]);
                        disp(' ')
                    end
                end
            else
                disp('No peaks calculated.')
            end
        end
        
        %% Print progress of Analysis
        function strCell=dispProgress(obj)
            strCell=cell(length(obj.progress),2);
            disp(' ')
            disp('Progression of Analysis:')
            for i=1:length(obj.progress)
                if ~isnan(obj.progress(i))
                    t=obj.progress(i);
                    tl='s';
                    if t>60
                        t=t/60;
                        tl='min';
                    end
                    if t>60
                        t=t/60;
                        tl='h';
                    end
                    
                    disp([obj.progresslabels{i} ' done in ' num2str(t,3) ' ' tl ])
                    strCell{i,1}=obj.progresslabels{i};
                    strCell{i,2}=[num2str(t,3) ' ' tl ];
                else
                    disp([obj.progresslabels{i} ' ...'])
                    strCell{i,1}=obj.progresslabels{i};
                    strCell{i,2}='...';
                end
            end
            disp(' ')
        end
        
        function obj=getSysInfo(obj)
            obj.sysInfo={};
            % CPU info
            if isunix
                [~,a]=system('cat /proc/cpuinfo | grep "model name" | uniq');
                obj.sysInfo{1}=a(1:end-1);
            else
                [~,a]=system('systeminfo');
                obj.sysInfo{1}=a;
            end
            
            % PoolSize info
            p = gcp('nocreate'); % If no pool, do not create new one.
            if isempty(p)
                obj.sysInfo{2} = 0;
            else
                obj.sysInfo{2} = p.NumWorkers;
            end
            
            %Version info
            obj.sysInfo{3}=obj.para.ver;
        end
        
        %% Load the exclude range into the gcanalysis object
        function obj=loadExcludeRanges(obj,exRange)
            try
                eval(['temp=' exRange ';']);
            catch
                error('loading exclude range: check format restrictions')
            end
            
            if isempty(temp)
                disp('exclude range is empty')
                obj.excludeRange=[];
                return
            end
            
            % check col number
            if size(temp,2)~=2
                error('loading exclude range: check format restrictions, more than two columns detected')
            end
            
            for i=1:size(temp,1)
                if temp(i,1)>=temp(i,2)
                    error(['loading exclude range: range #' num2str(i) 'is not increasing' ])
                end
            end
            
            disp('Exclude ranges:')
            for i=1:size(temp,1)
                disp(['  ' num2str(temp(i,1)) ' - ' num2str(temp(i,2)) ])
            end
            
            obj.excludeRange=temp;
        end
        
        %% Load the integration range into the gcanalysis object
        function obj=loadIntRanges(obj,intRange)
            try
                eval(['temp=' intRange ';']);
            catch
                error('loading integration range: check format restrictions')
            end
            
            if isempty(temp)
                disp('integration range is empty')
                obj.integrationRange=[];
                return
            end
            
            % check col number
            if size(temp,2)~=2
                error('loading integration range: check format restrictions, more than two columns detected')
            end
            
            for i=1:size(temp,1)
                if temp(i,1)>=temp(i,2)
                    error(['loading integration range: range #' num2str(i) 'is not increasing' ])
                end
            end
            
            disp('Integration ranges:')
            for i=1:size(temp,1)
                disp(['  ' num2str(temp(i,1)) ' - ' num2str(temp(i,2)) ])
            end
            
            obj.integrationRange=temp;
            
            %remove ranges outside of the time grid
            delidx=[];
            for i=1:size(obj.integrationRange,1)
                if obj.integrationRange(i,2)<min(obj.rawP.timegrid) || obj.integrationRange(i,1)>max(obj.rawP.timegrid)
                    delidx(end+1)=i;
                end
            end
            obj.integrationRange(delidx,:)=[];
        end
        
        
        %% Get exclude indices
        function [idxRangeBool]=getExcludeIndices(obj,gx)
            if nargin==1
                gx=obj.rawP.timegrid;
            end
            
            er=obj.excludeRange;
            
            idxRangeBool=true(size(gx));
            
            %no ranges present
            if isempty(er)
                return
            end
            
            for j=1:size(er,1)
                idxRangeBool(gx>=er(j,1) & gx<=er(j,2))=false;
            end
        end
        
        %% Get integration indices
        function [idxRangeBool]=getIntegrationIndices(obj,gx)
            if nargin==1
                gx=obj.rawP.timegrid;
            end
            
            er=obj.integrationRange;
            
            idxRangeBool=true(size(gx));
            
            %no ranges present
            if isempty(er)
                return
            end
            
            for j=1:size(er,1)
                idxRangeBool(gx>=er(j,1) & gx<=er(j,2))=false;
            end
        end
        
        %% Update cluster deadzones according to exclude ranges
        function updateClusterDeadzones(obj)
            
            er=obj.excludeRange;
            
            %no ranges present
            if isempty(er)
                return
            end
            
            for i=1:obj.rescl.anz_prob
                dz=obj.clusterP(i).deadzone;
                gx=obj.clusterP(i).timegrid;
                
                for j=1:size(er,1)
                    dz(gx>=er(j,1) & gx<=er(j,2))=false;
                end
                obj.clusterP(i).deadzone=dz;
            end
        end
        
        %% Calc integration ranges
        % the array of rangeInt objects is set here
        function calcIntegrationRanges(obj)
            if isempty(obj.integrationRange)
                return
            end
            
            % get global signal area
            area_global=trapz(obj.profPF.timegrid,obj.profPF.signal);
            
            obj.intRangeArray(size(obj.integrationRange,1),1)=rangeInt();
            
            for i=1:length(obj.intRangeArray)
                obj.intRangeArray(i).init(obj.integrationRange(i,:),...
                    obj.profPF.timegrid,...
                    obj.profPF.signal,...
                    area_global);
            end
            
        end
        
        %% Write peaks
        function writeIntRanges(obj,fpath)
            if isempty(obj.intRangeArray)
                return
            end
            
            if nargin==1
                fpath=[obj.para.tempf filesep 'int_res.csv'];
            end
            fid = fopen(fpath,'w');
            fprintf(fid,'%12s, ',...
                'index',...
                'max location',...
                'left bound',...
                'right bound',...
                'height',...
                'sig area');
            fprintf(fid,'%12s','rel sig area');
            fprintf(fid,'\n');
            for i=1:length(obj.intRangeArray)
                
                fprintf(fid,'%12f, ',i ,obj.intRangeArray(i).maxlocation, ...
                    obj.intRangeArray(i).lbound, ...
                    obj.intRangeArray(i).rbound, ...
                    obj.intRangeArray(i).height, ...
                    obj.intRangeArray(i).intinfo.sigAbs);
                fprintf(fid,'%12f',obj.intRangeArray(i).intinfo.sigRel);
                fprintf(fid,'\n');
            end
            fclose(fid);
        end
        
        
        
        
        %% Post processing peaks - Splitting
        function postprocSwapSubpeaks(obj)
            try
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Swap subpeaks - Start ###'])
                end
                
                % search for problemativ cases
                fullgx=obj.profPF.timegrid;
                ma_global=max(obj.profPF.signal);
                
                for i=1:length(obj.peaks)-1
                    %identify problematic cases
                    
                    %check if there is any overlap
                    rb1=obj.peaks(i).peakstruct.rbound;
                    lb2=obj.peaks(i+1).peakstruct.lbound;
                    if rb1<lb2
                        continue
                    end
                    
                    mx1=obj.peaks(i).peakstruct.maxlocation;
                    mx2=obj.peaks(i+1).peakstruct.maxlocation;
                    midx1=find(fullgx>=mx1,1,'first');
                    midx2=find(fullgx<=mx2,1,'last');
                    
                    %get tops of both peaks
                    mh1=obj.profPF.signal(midx1);
                    mh2=obj.profPF.signal(midx2);
                    
                    %check if the peaks are large enough
                    if min(mh1,mh2)<10*obj.rescl.noiselvl
                        continue
                    end
                    
                    %ensure that there's at least some intermediate point
                    if midx1+5<midx2
                        
                        %get minimum location between tops
                        interGrid=obj.profPF.timegrid(midx1:midx2);
                        interSig=obj.profPF.signal(midx1:midx2);
                        [minH,minIdx]=min(interSig);
                        minPos=interGrid(minIdx);
                        
                        % check if there is a sufficiently deep valley between
                        % peak maxima
                        if minH > 0.8*min(mh1,mh2)
                            continue
                        end
                        
                    else
                        continue
                    end
                    
                    %get swap peak candidates
                    swap1=false(1,obj.peaks(i).modelnumber);
                    swap2=false(1,obj.peaks(i+1).modelnumber);
                    
                    % left peak
                    for j=1:obj.peaks(i).modelnumber
                        if obj.peaks(i).parameters(1,j)>minPos
                            swap1(j)=true;
                        end
                    end
                    
                    % right peak
                    for j=1:obj.peaks(i+1).modelnumber
                        if obj.peaks(i+1).parameters(1,j)<minPos
                            swap2(j)=true;
                        end
                    end
                    
                    
                    if max(swap1) || max(swap2)
                        if obj.para.info>=2
                            disp(['=>  Peak swap for peaks ' num2str(i) ' and ' num2str(i+1)])
                        end
                        
                        % swap peaks between peak.parameters matrices
                        %get old parameter matrices
                        paranew1=obj.peaks(i).parameters;
                        paranew2=obj.peaks(i+1).parameters;
                        %remove peaks that should be swapped
                        paranew1(:,swap1)=[];
                        paranew2(:,swap2)=[];
                        %add the swap peaks to the other parameter matrix
                        paranew1=[paranew1 obj.peaks(i+1).parameters(:,swap2)];
                        paranew2=[obj.peaks(i).parameters(:,swap1) paranew2];
                        
                        % restore order by subpeak center location
                        [~,order]=sort(paranew1(1,:));
                        paranew1=paranew1(:,order);
                        [~,order]=sort(paranew2(1,:));
                        paranew2=paranew2(:,order);
                        
                        %assign new parameters
                        obj.peaks(i).parameters=paranew1;
                        obj.peaks(i).modelnumber=size(paranew1,2);
                        obj.peaks(i+1).parameters=paranew2;
                        obj.peaks(i+1).modelnumber=size(paranew2,2);
                        
                        % reevaluate peak models (no correction present)
                        obj.peaks(i).updatePropAfterCorr(fullgx,ma_global);
                        obj.peaks(i+1).updatePropAfterCorr(fullgx,ma_global);
                        
                        
                    end
                end
                
                
                if obj.para.info>=1
                    disp(['### Postprocess peaks - Swap subpeaks - End ###'])
                end
                
            catch ex
                obj.errorLog(ex,'postprocSwapSubpeaks')
            end
        end
        
        %% Post processing peaks - Remove Outlier Subpeaks
        function removeOutlierSubpeaks(obj)
            %remove subpeaks that are located outside of the range between
            %lbound and rbound
            
            try
                fullgx=obj.profPF.timegrid;
                ma_global=max(obj.profPF.signal);
                
                for i=1:length(obj.peaks)
                    del=false(1,obj.peaks(i).modelnumber);
                    
                    %identify outlier subpeaks
                    lb=obj.peaks(i).peakstruct.lbound;
                    rb=obj.peaks(i).peakstruct.rbound;
                    for j=1:obj.peaks(i).modelnumber
                        if obj.peaks(i).parameters(1,j) < lb || obj.peaks(i).parameters(1,j) > rb
                            del(j)=true;
                        end
                    end
                    
                    %remove outlier
                    obj.peaks(i).parameters(:,del)=[];
                    obj.peaks(i).modelnumber=size(obj.peaks(i).parameters,2);
                    
                    % reevaluate peak model
                    obj.peaks(i).updatePropAfterCorr(fullgx,ma_global);
                end
            catch ex
                obj.errorLog(ex,'removeOutlierSubpeaks')
            end
        end
        
        %% Create errorLogFile and save object
        function errorLog(obj,ex,comment)
            %ERRORLOG saves error and current status
            
            save_file = [obj.para.tempf 'ErrorStatus.mat'];
            save(save_file)
            
            fid = fopen(fullfile(obj.para.tempf, 'ErrorLogFile.txt'), 'w');
            fprintf(fid, '%s: %s\n','File', obj.para.tempf);
            fprintf(fid, '%s: %s\n','Ver.', obj.para.ver);
            if nargin >= 3 
                fprintf(fid, '%s\n', comment);
            end
            
            if nargin >= 2
                fprintf(fid, '%s\n', ex.getReport);
            else
                fprintf(fid, 'Error Report by user.\n');
            end
            fprintf(fid, '%s: %s\n', 'current status saved as', save_file);
            
            %error(sprintf('very long error msg\nbla blah bla'))
            error(sprintf(['Error-information can be found in the results-folder:\n',...
                [save_file '\n'],...
                'You may contact the developers under:\n',...
                'martina.beese@uni-rostock.de\n',...
                'tomass.andersons@uni-rostock.de\n',...
                'Your request should contain the ErrorLogFile.txt and the ErrorStatus.mat ']))
            
        end
        
    end
    %method end
    
    
end
%class end








































