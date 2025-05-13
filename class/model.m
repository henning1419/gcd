classdef model < profileC
    % set of peaks used in the optimization process
    % obtained by clusters from GCanalysis.calcCluster;
    
    properties
        % timegrid - inherited from profile - double vector
        % time subgrid for peak cluster
        
        % signal - inherited from profile -double vector
        % signal for peak cluster
        
        % set of peaks(class)
        peaks=[];
        
        % temporary model during optimization
        modeleval
        
        %scal (DEL???)
        
        % deadzone - bool vector(length(timegrid))
        % indicates in which regions on the timegrid peak fitting is allowed.
        % regions are deactivated, if e.g. an optimization fails or the
        % fitting error is to high.
        % (true = allowed, false = forbidden)
        deadzone
        
        % noCurveAddCounter - scalar integer
        % counts failed peak fits, too many will trigger the deactivation
        % of the cluster in the optimization
        noCurveAddCounter=0;
        
        % resfi - struct - relic from early version - generated in GCanalysis.calcFitPar
        %   peaks : same as model.peaks
        %    scal : 1
        %   zuord : assignment of overlapping subpeaks to "real" peaks - calculated by calcZuordAlt2 
        resfi=[];
    end
    
    methods
        %% CONSTRUCTOR
        function obj = model(timegrid,signal)
            if nargin==2
                obj.timegrid = timegrid;
                obj.signal = signal;

                obj.modeleval=zeros(size(signal));
                obj.deadzone=true(size(signal));
                
            elseif nargin==0
                
            end
        end
        
        %% DEPRECATED
%         function resfi=calcModel(obj,para,noiselvl)
%             d=para.d;
%             M=para.M;
%             N=2*M+1;
%             
%             resfi=corr_fit2(para,obj.timegrid,obj.signal,d,N,noiselvl);
%             obj.modeleval=resfi.gcmod*resfi.scal;
%             
%         end
        
        %% Post fit - combination of single sub peaks to combined peaks (alternative)
        function Zuord=calcZuordAlt2(obj,para)
         % Assign small subpeaks to the larger ones
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % is the intersecting area of a curve larger as p_corell*100 percent
            % than it is assigned to the larger one.
            % modification: subpeaks are always assigned to those subpeaks
            % with the largest overlap area
            
            if para.info>=2
                disp(' ')
                disp('=> Assignement')
            end
            
            %get all sub peak parameters columnwise
            Pm=zeros(6,length(obj.peaks));
            for i=1:length(obj.peaks)
                Pm(:,i)=obj.peaks{i}.para;
            end
            
            % init the assignement structure:
            % in the beginning Zuord{i} contains only the i-th subpeak
            % in the end Zuord{i} may:
            % 1) be emtpy - then the subpeak has been assigned to another subpeak
            % 2) contain one or more suboeak indices. together they form
            %    the combined peak
            
            Zuord=cell(length(obj.peaks),1);
            for i=1:length(obj.peaks)
                Zuord{i}=i;
            end
            
            
            % make assignements
            ii=1;
            doneAssign=1;
            while doneAssign
                doneAssign=0;
                if para.info>=2
                    disp(['<> Overlap run ' num2str(ii)])
                    ii=ii+1;
                end
                
                %sort by peak height
                peakma=zeros(1,length(obj.peaks));
                for i=1:length(obj.peaks)
                    peakma(i)=obj.peaks{i}.para(4);
                end
                [~,sorting]=sort(peakma);
                
                
                for i=sorting
                    % get current sub peak info
                    if ~isempty(Zuord{i})
                        cZu=Zuord{i};
                        mod=mod_ga(obj.timegrid,Pm(:,cZu(1)));
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
                            mod=mod+mod_ga(obj.timegrid,Pm(:,cZu(k)));
                        end
                    else
                        continue
                    end
                    
                    %compare current peak with all remaining
                    
                    %vector of overlap
                    vOver=nan(1,length(obj.peaks));
                    vAssign=cell(1,length(obj.peaks));
                    for j=1:length(obj.peaks)
                        if i~=j
                            %get sub peak info
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
                                
                                mod2=mod_ga(obj.timegrid,Pm(:,cZu2(1)));
                                for k=2:length(cZu2)
                                    mod2=mod2+mod_ga(obj.timegrid,Pm(:,cZu2(k)));
                                end
                                
                                I_t=trapz(obj.timegrid,mod);
                                I_t2=trapz(obj.timegrid,mod2);
                                I_inter=trapz(obj.timegrid,min(mod,mod2));
                                [m1,i1]=max(mod);
                                [m2,i2]=max(mod2);
                                modc=mod+mod2;
                                mc1=modc(i1);
                                mc2=modc(i2);
                                min_mc=min([mc1,mc2]);
                                mic=min(modc(min(i1,i2):max(i1,i2)));
                                
                                
                                %Peaks auf Tailing verhindern
                                %Wenn Tal zwischen zwei Maxima tief genug
                                %dann ueberspringe
                                
                                %cancel criteria
                                if abs((min_mc-mic)/min_mc)>=0.0045 && ...
                                        abs((I_t-I_inter)./I_inter)>0.005 && ...
                                        abs((I_t2-I_inter)./I_inter)>0.005
                                    continue
                                end
                                
                                %accept criteria
                                if I_inter/I_t>para.peakoverlap || I_inter/I_t2>para.peakoverlap
                                    %check if there is another peak with a
                                    %larger peak overlap
                                    
                                    vAssign{j}=Zuord{j};
                                    vOver(j)=max(I_inter/I_t,I_inter/I_t2);
%                                     Zuord{i}=[Zuord{i} Zuord{j}];
%                                     Zuord{j}=[];
                                    doneAssign=1;
                                    continue
                                    
                                end
                                
                                %only check if not already accepted
                                if i1<i2 && m1>m2*4 && I_inter/I_t2>0.1*para.peakoverlap && abs((mc2-mic)/mc2)<0.02 
                                    if para.info>1
                                        disp('<> Tail found/appended')
                                    end
                                    vAssign{j}=Zuord{j};
                                    vOver(j)=max(I_inter/I_t,I_inter/I_t2);
%                                     Zuord{i}=[Zuord{i} Zuord{j}];
%                                     Zuord{j}=[];
                                    doneAssign=1;
                                end
                                
                            end
                        end
                    end
                    
                    % overlapping peaks have been found
                    % assign those with highest overlap
                    if doneAssign
                        
                        [maVOver,tempidx]=max(vOver);
                        if isnan(maVOver)
                            continue;
                        end
                        Zuord{i}=[Zuord{i} Zuord{tempidx}];
                        Zuord{tempidx}=[];
                    end
                    
                end
            end
        end
        
        %% Post fit - combination of single sub peaks to combined peaks (alternative)
        function Zuord=calcZuordAlt(obj,para)
            % Assign small subpeaks to the larger ones
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % is the intersecting area of a curve larger as p_corell*100 percent
            % than it is assigned to the larger one.
            % modification: subpeaks are always assigned to those subpeaks
            % with the largest overlap area
            
            if para.info>=2
                disp(' ')
                disp('=> Assignement')
            end
            
            %get all sub peak parameters columnwise
            Pm=zeros(6,length(obj.peaks));
            for i=1:length(obj.peaks)
                Pm(:,i)=obj.peaks{i}.para;
            end
            
            % init the assignement structure:
            % in the beginning Zuord{i} contains only the i-th subpeak
            % in the end Zuord{i} may:
            % 1) be emtpy - then the subpeak has been assigned to another subpeak
            % 2) contain one or more suboeak indices. together they form
            %    the combined peak
            
            Zuord=cell(length(obj.peaks),1);
            for i=1:length(obj.peaks)
                Zuord{i}=i;
            end
            
            
            % make assignements
            ii=1;
            doneAssign=1;
            while doneAssign
                doneAssign=0;
                if para.info>=2
                    disp(['<> Overlap run ' num2str(ii) ' of ' num2str(para.overlap_runs)])
                    ii=ii+1;
                end
                
                for i=1:length(obj.peaks)
                    % get current sub peak info
                    if ~isempty(Zuord{i})
                        cZu=Zuord{i};
                        mod=mod_ga(obj.timegrid,Pm(:,cZu(1)));
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
                            mod=mod+mod_ga(obj.timegrid,Pm(:,cZu(k)));
                        end
                    else
                        continue
                    end
                    
                    %compare current peak with all remaining
                    
                    %vector of overlap
%                     vOver=nan(1,length(obj.peaks));
%                     vAssign=cell(1,length(obj.peaks));
                    for j=1:length(obj.peaks)
                        if i~=j
                            %get sub peak info
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
                                
                                mod2=mod_ga(obj.timegrid,Pm(:,cZu2(1)));
                                for k=2:length(cZu2)
                                    mod2=mod2+mod_ga(obj.timegrid,Pm(:,cZu2(k)));
                                end
                                
                                I_t=trapz(obj.timegrid,mod);
                                I_t2=trapz(obj.timegrid,mod2);
                                I_inter=trapz(obj.timegrid,min(mod,mod2));
                                [m1,i1]=max(mod);
                                [m2,i2]=max(mod2);
                                modc=mod+mod2;
                                mc1=modc(i1);
                                mc2=modc(i2);
                                min_mc=min([mc1,mc2]);
                                mic=min(modc(min(i1,i2):max(i1,i2)));
                                
                                
                                %Peaks auf Tailing verhindern
                                %Wenn Tal zwischen zwei Maxima tief genug
                                %dann ueberspringe
                                
                                %cancel criteria
                                if abs((min_mc-mic)/min_mc)>=0.0045 && ...
                                        abs((I_t-I_inter)./I_inter)>0.005 && ...
                                        abs((I_t2-I_inter)./I_inter)>0.005
                                    continue
                                end

                                %accept criteria
                                if I_inter/I_t>para.peakoverlap || I_inter/I_t2>para.peakoverlap
                                    %check if there is another peak with a
                                    %larger peak overlap
                                    
                                    %%
                                    okValiTest=1;
                                    for jjj=1:length(obj.peaks)
                                        if j~=jjj && jjj~=i
                                            %get sub peak info
                                            if ~isempty(Zuord{jjj})
                                                cZuVali=Zuord{jjj};
                                                
                                                cce=Pm(1,cZuVali(1));
                                                
                                                csiL=Pm(2,cZuVali(1));
                                                csiR=Pm(3,cZuVali(1));
                                                cleftb=cce-3*csiL;
                                                crightb=cce+3*csiR;
                                                for k=2:length(cZuVali)
                                                    cce=Pm(1,cZuVali(k));
                                                    csiL=Pm(2,cZuVali(k));
                                                    csiR=Pm(3,cZuVali(k));
                                                    cleftb=min(cleftb,  cce-3*csiL );
                                                    crightb=max(crightb,  cce+3*csiR );
                                                end
                                                
                                                %fast check, trapz, interp1 etc cost time
                                                if leftb>crightb || cleftb >rightb
                                                    continue
                                                end
                                                
                                                modVali=mod_ga(obj.timegrid,Pm(:,cZuVali(1)));
                                                for k=2:length(cZuVali)
                                                    modVali=modVali+mod_ga(obj.timegrid,Pm(:,cZuVali(k)));
                                                end
                                                
                                                %compare I_t2 with I_Vali
                                                I_Vali=trapz(obj.timegrid,modVali);
                                                I_interVali=trapz(obj.timegrid,min(mod2,modVali));
                                                
                                                if max(I_interVali/I_t2,I_interVali/I_Vali) > 1.01*max(I_inter/I_t,I_inter/I_t2)
                                                    okValiTest=0;
                                                    if 0
                                                        figure(99);
                                                        clf
                                                        plot(obj.timegrid,mod)
                                                        hold on
                                                        plot(obj.timegrid,mod2)
                                                        plot(obj.timegrid,modVali)
                                                        legend('current','test','cali')
                                                        1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                    %%
                                    
                                    if okValiTest
%                                         vAssign{j}=Zuord{j};
%                                         vOver(j)=max(I_inter/I_t,I_inter/I_t2);
                                        Zuord{i}=[Zuord{i} Zuord{j}];
                                        Zuord{j}=[];
                                        doneAssign=1;
                                        continue
                                    end
                                    
                                end
                                
                                %only check if not already accepted
                                if i1<i2 && m1>m2*4 && I_inter/I_t2>0.1*para.peakoverlap && abs((mc2-mic)/mc2)<0.02 
                                    if para.info>1
                                        disp('<> Tail found/appended')
                                    end
%                                     vAssign{j}=Zuord{j};
%                                     vOver(j)=max(I_inter/I_t,I_inter/I_t2);
                                    Zuord{i}=[Zuord{i} Zuord{j}];
                                    Zuord{j}=[];
                                    doneAssign=1;
                                end
                                
                            end
                        end
                    end
                    
                    % overlapping peaks have been found
                    % assign those with highest overlap
%                     if doneAssign
%                         %check if 
%                         
%                         [~,tempidx]=max(vOver);
%                         Zuord{i}=[Zuord{i} Zuord{tempidx}];
%                         Zuord{midx}=[];
%                     end
                    
                end
            end
        end
        
        
        %% Post fit - combination of single sub peaks to combined peaks
        function Zuord=calcZuord(obj,para)
            
            % Assign small subpeaks to the larger ones
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % is the intersecting area of a curve larger as p_corell*100 percent
            % than it is assigned to the larger one.
            
            if para.info>=2
                disp(' ')
                disp('=> Assignement')
            end
            
            %get all sub peak parameters columnwise
            Pm=zeros(6,length(obj.peaks));
            for i=1:length(obj.peaks)
                Pm(:,i)=obj.peaks{i}.para;
            end
            
            % init the assignement structure:
            % in the beginning Zuord{i} contains only the i-th subpeak
            % in the end Zuord{i} may:
            % 1) be emtpy - then the subpeak has been assigned to another subpeak
            % 2) contain one or more suboeak indices. together they form
            %    the combined peak
            
            Zuord=cell(length(obj.peaks),1);
            for i=1:length(obj.peaks)
                Zuord{i}=i;
            end
            
            
            % make assignements
            for ii=1:para.overlap_runs
                if para.info>=2
                    disp(['<> Overlap run ' num2str(ii) ' of ' num2str(para.overlap_runs)])
                end
                
                for i=1:length(obj.peaks)
                    % get current sub peak info
                    if ~isempty(Zuord{i})
                        cZu=Zuord{i};
                        mod=mod_ga(obj.timegrid,Pm(:,cZu(1)));
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
                            mod=mod+mod_ga(obj.timegrid,Pm(:,cZu(k)));
                        end
                    else
                        continue
                    end
                    %compare current peak with all remaining
                    for j=1:length(obj.peaks)
                        if i~=j
                            %get sub peak info
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
                                
                                mod2=mod_ga(obj.timegrid,Pm(:,cZu2(1)));
                                for k=2:length(cZu2)
                                    mod2=mod2+mod_ga(obj.timegrid,Pm(:,cZu2(k)));
                                end
                                %                 mod2=mod_ga(gx,Pm(:,j));
                                I_t=trapz(obj.timegrid,mod);
                                I_t2=trapz(obj.timegrid,mod2);
                                I_inter=trapz(obj.timegrid,min(mod,mod2));
                                [m1,i1]=max(mod);
                                [m2,i2]=max(mod2);
                                modc=mod+mod2;
                                mc1=modc(i1);
                                mc2=modc(i2);
                                min_mc=min([mc1,mc2]);
                                mic=min(modc(min(i1,i2):max(i1,i2)));
                                
                                
                                %Peaks auf Tailing verhindern
                                %Wenn Tal zwischen zwei Maxima tief genug
                                %dann ueberspringe
                                
                                %cancel criteria
                                if abs((min_mc-mic)/min_mc)>=0.0045 && ...
                                        abs((I_t-I_inter)./I_inter)>0.005 && ...
                                        abs((I_t2-I_inter)./I_inter)>0.005
                                    continue
                                end
                                
                                %accept criteria
                                if I_inter/I_t>para.peakoverlap || I_inter/I_t2>para.peakoverlap
                                    Zuord{i}=[Zuord{i} Zuord{j}];
                                    Zuord{j}=[];
                                    continue
                                end
                                
                                
                                if i1<i2 && m1>m2*4 && I_inter/I_t2>0.1*para.peakoverlap && abs((mc2-mic)/mc2)<0.02
                                    if para.info>1
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
        end
        
        
        
        %%
        function b=checkClusterAct(obj,para,id,noiselvl)
            b=1;
            % too many peak fit tries without success
            if obj.noCurveAddCounter >=para.retry
                b=0;
                if para.info>=1
                    disp(['<> Set Cluster ' num2str(id) ' to inactive: too many failed fits'])
                end
                return
            end
            
            % deadzone spans complete signal
            if max(obj.deadzone)==0
                b=0;
                if para.info>=1
                    disp(['<> Set Cluster ' num2str(id) ' to inactive: deadzone spans complete signal'])
                end
                return
            end
            
            % reached maximal number of modelled peaks
            if length(obj.peaks)>=para.maxpeak_nr
                b=0;
                if para.info>=1
                    disp(['<> Set Cluster ' num2str(id) ' to inactive: max peak number reached'])
                end
                return
            end
            
            % fit is well enough
            res=obj.signal-obj.modeleval;
            err=abs(trapz(res))/abs(trapz(obj.signal));
            if err<=para.errTol
                b=0;
                if para.info>=1
                    disp(['<> Set Cluster ' num2str(id) ' to inactive: fit error tolerance reached'])
                end
                return
            end
            
            % largest peak in residue is below x*noiselvl
            if max(res)<=3*noiselvl
                b=0;
                if para.info>=1
                    disp(['<> Set Cluster ' num2str(id) ' to inactive: max(residue) is below tolerance'])
                end
                return
            end
            
        end
        
        function [idx,smooth_sig_max]=getNextPeakCalc(obj,currCalcZone,para)
            d=para.d;
            M=para.M;
            N=2*M+1;
            
            signal=obj.signal;
            mod=obj.modeleval;
            
            % no free peak fit area available
            if max(~currCalcZone & obj.deadzone)==0
                idx=-1;
                if para.info>=2
                    disp(['<> getNextPeakCalc warning: no free peak fit area available'])
                end
                smooth_sig_max=[];
                return
            end
            
            currGC=signal-mod;
            [d_smooth,~,~,~]= sgfilt(d,N,currGC);
            %get idx for max model error
            [smooth_sig_max,idxred]=max(d_smooth(~currCalcZone & obj.deadzone));

            %get full range of indices
            vfull=1:length(signal);
            %reduce by cutting irrelevant ranges
            vreduce=vfull(~currCalcZone & obj.deadzone);
            %obtain idx on full timegrid/signal
            idx=vreduce(idxred);
            
            debug=0;
            if debug
                figure(1);
                clf
                subplot(4,1,1);
                plot(obj.timegrid,obj.signal);
                hold on
                plot(obj.timegrid,obj.modeleval,'r--','linewidth',2);
                plot(obj.timegrid(idx),obj.modeleval(idx),'ro')
                title('sig + mod')
                
                subplot(4,1,2)
                plot(obj.timegrid,currGC);
                hold on
                tg=obj.timegrid(~currCalcZone);
                plot(tg, zeros(size(tg)),'bx')
                plot(obj.timegrid(idx),obj.modeleval(idx),'ro')
                title('calc zone')
                
                subplot(4,1,3)
                plot(obj.timegrid,currGC);
                hold on
                tg=obj.timegrid(obj.deadzone);
                plot(tg, zeros(size(tg)),'rx')
                plot(obj.timegrid(idx),obj.modeleval(idx),'ro')
                title('dead zone')
                
                subplot(4,1,4)
                plot(obj.timegrid,currGC);
                hold on
                tg=obj.timegrid(~currCalcZone & obj.deadzone);
                plot(tg, zeros(size(tg)),'gx')
                plot(obj.timegrid(idx),obj.modeleval(idx),'ro')
                title('calc and dead zone')
                1;
            end
        end
        
        function p=getApproxPeak(obj,timegridIdx,smooth_sig_max,para)
            p.success=0;
            signal=obj.signal;
            modev=obj.modeleval;
            resev=signal-modev;
            timegrid=obj.timegrid;
            info=para.info;
            
            %tg_idx=timegrid(timegridIdx); %old Vxmin
            %sig_idx=smooth_sig_max; %old
            
            % get local maximum in signal
            try
                % search left side
                sig_loc_maxL=smooth_sig_max;
                jL=timegridIdx;
                if jL>=2
                    while sig_loc_maxL<resev(jL-1)
                        sig_loc_maxL=resev(jL-1);
                        jL=jL-1;
                        if jL==1
                            break
                        end
                    end
                end
                
                % search right side
                sig_loc_maxR=smooth_sig_max;
                jR=timegridIdx;
                if jR<=length(resev)-1
                    while sig_loc_maxR<resev(jR+1)
                        sig_loc_maxR=resev(jR+1);
                        jR=jR+1;
                        if jR==length(resev)
                            break
                        end
                    end
                end
                
                % combine results of left and right search
                if sig_loc_maxL>=sig_loc_maxR
                    centerIdx=jL;
                    height=sig_loc_maxL;
                else
                    centerIdx=jR;
                    height=sig_loc_maxR;
                end
                center=timegrid(centerIdx);
            catch
                disp('<> getApproxPeak error: ID 0')
            end
            
            % get width by approximate when the peak reaches 3/5 of iths height
            idxL=find(((timegrid<=center)+(   (resev)<=height*3/5   ))==2,1,'last'); %3/5
            widthL=center-timegrid(idxL);
            idxR=find(((timegrid>=center)+(   (resev)<=height*3/5   ))==2,1,'first'); %3/5
            widthR=timegrid(idxR)-center;
            
            if isempty(widthL) || isempty(widthR)
                widthid=1;
                if info>1
                    disp(['<> getApproxPeak: widthL/R empty'])
                end
            else
                if widthL<widthR
                    xid=idxL;
                    width=widthL;
                else
                    xid=idxR;
                    width=widthR;
                end
                widthid=abs(centerIdx-xid);
            end
            % ignore to narrow peaks, extend cluster deadzone
            if widthid<=2
                if info>1
                    disp(['<> getApproxPeak continue due to too narrow peak: widthid=' num2str(widthid)])
                end
                obj.deadzone(max(timegridIdx-ceil(para.bw*2),1):min(timegridIdx+floor(para.bw*2),length(signal))  )=false;
                return
            end
            
            % tight peak range (constraint model=signal)
            le_idx=max([1 centerIdx-floor(0.5*widthid)]);
            ri_idx=min([length(resev) centerIdx+ceil(0.5*widthid)]);
            
            % far peak range (constraint model<=signal)
            le_idx_w=max([1 centerIdx-10*widthid]);
            ri_idx_w=min([length(resev) centerIdx+10*widthid]);
            
            % get active signal ranges for constraint evaluation
            % dif -> model<=signal; neg -> model=signal
            act_dif=false(1,length(resev));
            act_dif(le_idx:ri_idx)=1;
            act_neg=false(1,length(resev));
            act_neg(le_idx_w:ri_idx_w)=1;
            
            idx_difStartInNeg=le_idx-le_idx_w+1;
            idx_difLen=ri_idx-le_idx+1;
            
            
            p.center=center;
            p.centerIdx=centerIdx;
            p.width=width;
            p.height=height;
            p.act_dif=act_dif;
            p.act_neg=act_neg;
            p.idx_difStartInNeg=idx_difStartInNeg;
            p.idx_difLen=idx_difLen;
            p.success=1;
            
            debug=0;
            if debug
                disp('## par')
                disp(['act_dif: lidx ' num2str(le_idx) ' ridx ' num2str(ri_idx) ])
                disp(['act_neg: lidx ' num2str(le_idx_w) ' ridx ' num2str(ri_idx_w) ])
                
                figure(5)
                clf
%                 subplot(2,1,1)
                plot(timegrid,resev);
                hold on
                plot(center,resev(centerIdx),'xr')
                plot(timegrid(idxL),resev(idxL),'>r')
                plot(timegrid(idxR),resev(idxR),'<r')
                title(['par ' num2str(centerIdx) ' of ' num2str(length(resev))])
                1;
            end
        end
        
        
        function j=prepareJob(obj,p,noiselvl)
            %shift time grid, such that center -> 0
            timegrid_shift=obj.timegrid-p.center;
            j.timegrid_shi_red=timegrid_shift(p.act_neg);
            
            
            %initial parameter
            p0=[0,...              center
                1.3*p.width,...      sigma left  (1.3*)
                1.3*p.width,...      sigma right (1.3*)
                0.95,...           height
                0.3, ...           linearity left
                0.3, ...           linearity right
                ];
            
            %lower/upper bounds
            lb=zeros(size(p0));
            ub=inf(size(p0));
            lb(1)=min(timegrid_shift(p.act_dif));
            ub(1)=max(timegrid_shift(p.act_dif));
            lb(2)=0.05*p0(2);
            lb(3)=0.05*p0(3);
            ub(2)=2.5*p0(2);
            ub(3)=10*p0(3);
            ub(5)=1;
            ub(6)=1;
            lb=lb(:);
            ub=ub(:);
            p0=p0(:);
            
            j.p0=p0;
            j.ub=ub;
            j.lb=lb;
            
            %get remaining signal
            res=obj.signal-obj.modeleval;
            res_reduce=res(p.act_neg);
            
            %scal to max 1 (cf. p0(4) = height)
            j.res_red_scal=max(res_reduce);
            j.res_red_scaled=res_reduce./j.res_red_scal;
            
            % get calczone from noiselvl -> peak -> noiselvl
            % 1 = calc, 0 = no calc
            idCZ_L=find(timegrid_shift<=0 & res<=2*noiselvl,1,'last');
            idCZ_R=find(timegrid_shift>=0 & res<=2*noiselvl,1,'first');
            if isempty(idCZ_L)
                idCZ_L=1;
            end
            if isempty(idCZ_R)
                idCZ_R=length(timegrid_shift);
            end
            j.calcZone=false(size(timegrid_shift));
            j.calcZone(idCZ_L:idCZ_R)=true;
            
            debug=0;
            if debug
                figure(5)
                xxl=xlim;
                plot(xxl,noiselvl*[1,1],'r--');
                1;
            end
        end
        
        
        function plot(obj,col)
            plot(obj.timegrid,obj.signal,'k')
            hold on
            if nargin==1
                plot(obj.timegrid,obj.modeleval,'b')
            elseif nargin==2
                plot(obj.timegrid,obj.modeleval,col)
            end
        end
    end
end

