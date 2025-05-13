classdef peak < profileC
    %PEAK Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % timegrid - inherited from profile - double vector
        % time subgrid between peakstruct.lbound and peakstruct.rbound
        
        % signal - inherited from profile -double vector
        % signal between peakstruct.lbound and peakstruct.rbound
        
        % modelnumber - integer
        % number of used subpeaks 
        modelnumber
        
        % parameters - double matrix(6 x modelnumber)
        % each column contains a parameter vector for a subpeak
        % entry - meaning
        %     1 - center (>=0)
        %     2 - left sigma (>=0)
        %     3 - right sigma (>=0)
        %     4 - height (>=0)
        %     5 - left linearity coefficient ([0,1]), 0=gauss, 1=linear(with smoothed ends) 
        %     6 - right linearity coefficient ([0,1])
        parameters
        
        % peakstruct - struct with fields 
        %     peakX_l2r: (deprecated)
        %  fullpeak_l2r: signal between lbound and rbound
        %       paramat: (deprecated) use peak.parameters instead
        %        lbound: left border of peak (first above peakMax/300)
        %        rbound: right border of peak (last above peakMax/300)
        %   maxlocation: location of the peak maximum
        %        height: peak height
        %     relheight: relative peak height compared to signal maximum height
        %   modelnumber: (deprecated) use peak.modelnumber instead
        peakstruct
        
        % correction - double vector
        % model correction between peakstruct.lbound and peakstruct.rbound
        correction
        
        % overlap/fit - double scalar
        % indicators for overlapping with other peaks
        % updated by GCanalysis.evalPeaks;
        overlap=0;
        fit=0;
        
        %intinfo - struct - contains integral information of the peak
        % -> intinfo.abs        - absolute model integrals
        % -> intinfo.rel        - relative model integrals
        % -> intinfo.sigAbs     - absolute signal based integrals
        % -> intinfo.sigRel     - relative signal based integrals
        % -> intinfo.agiAbs     - absolute agilent based integrals
        % -> intinfo.agiRel     - relative agilent based integrals
        intinfo=[];
        
    end
    
    methods
        function obj = peak(timegrid,signal,modelnumber,parameters,peakstruct)
            if nargin==5
                obj.timegrid=timegrid;
                obj.signal=signal;
                obj.modelnumber=modelnumber;
                obj.parameters=parameters;
                obj.peakstruct=peakstruct;
                obj.correction=zeros(size(timegrid));
            elseif nargin==0
                obj.timegrid=[];
                obj.signal=[];
                obj.modelnumber=0;
                obj.parameters=[];
                obj.peakstruct=[];
                obj.correction=[];
            end
        end
        
        function b=istailing(obj)
            lb=obj.peakstruct.lbound;
            ce=obj.peakstruct.maxlocation;
            rb=obj.peakstruct.rbound;
            hi=obj.peakstruct.height;
            if obj.modelnumber==1
                b=0;
                return
            end
            b=0;
            
            %ration of left flank to right flank
            if (rb-ce)/(ce-lb)>3
                b=1;
            end
        end
        
        function b=istailingDetail(obj,sigP,noiselvl)
            lb=obj.peakstruct.lbound;
            ce=obj.peakstruct.maxlocation;
            rb=obj.peakstruct.rbound;
            hi=obj.peakstruct.height;
            if obj.modelnumber==1
                b=0;
                return
            end
            
            
            gx=sigP.timegrid;
            gc=sigP.signal;
            
            b=0;
            %ration of left flank to right flank
            if (rb-ce)/(ce-lb)>3 && hi>=40*noiselvl
                idx_ce=find(gx>=ce,1,'first');
                if gc(idx_ce) <= 1.05*hi
                    b=1;
                end
            end
        end
        
        %%
        function mod=getmod(obj,gx)
            mod=zeros(length(gx),1);
            param=obj.parameters;
            for jj=1:size(param,2)
                mod=mod+mod_ga(gx,param(:,jj));
            end
            if ~isempty(obj.correction)
                mod=mod+interp1(obj.timegrid,obj.correction,gx,'linear',0);
            end
            
        end
        
        %%
        function [pArray]=splitMulti(obj,spxArray,gx,ma_global,info)
            if nargin<6
                info=1;
            end
            
            if isempty(spxArray)
                error('splitting array is empty')
            end
            
            spxArray=sort(spxArray);
            pArray=peak.empty;
            
            [pL,pR]=obj.split(spxArray(1),gx,ma_global,info);
            pArray(end+1)=pL;
            
            for i=2:length(spxArray)
                [pL,pR]=pR.split(spxArray(2),gx,ma_global,info);
                pArray(end+1)=pL;
            end
            pArray(end+1)=pR;
            
        end
        
        %%
        function [pl,pr]=split(obj,spx,gx,ma_global,info)
            if nargin<6
                info=1;
            end
            %TODO correction splitten
            paral=obj.parameters(:,obj.parameters(1,:)<=spx);
            parar=obj.parameters(:,obj.parameters(1,:)>spx);
            modelnl=size(paral,2);
            modelnr=size(parar,2);
            
            %Left
            fullpeak=zeros(length(gx),1);
            for i=1:modelnl
                fullpeak=fullpeak+mod_ga(gx,paral(:,i));
            end
            
            [ma,id]=max(fullpeak);
            left=find(fullpeak>ma/300,1,'first');
            right=find(fullpeak>ma/300,1,'last');
            if info
                disp('Split Peak Left        ');
                disp(['    Left bound:        ' num2str(gx(left))]);
                disp(['    Peak max location: ' num2str(gx(id))]);
                disp(['    Right bound:       ' num2str(gx(right))]);
                disp(['    Peak height:       ' num2str(ma)]);
                disp(['    Peak rel. height:  ' num2str(ma/ma_global)]);
                disp(['    Model number:      ' num2str(modelnl)]);
            end
            p.peakX_l2r=gx(left:right);
            p.fullpeak_l2r=fullpeak(left:right);
            p.paramat=paral;
            p.lbound=gx(left);
            p.rbound=gx(right);
            p.maxlocation=gx(id);
            p.height=ma;
            p.relheight=ma/ma_global;
            p.modelnumber=modelnl;
            
            pl=peak(p.peakX_l2r,p.fullpeak_l2r,modelnl,paral,p);
            %             pl.correction=obj.correction()
            
            %Right
            fullpeak=zeros(length(gx),1);
            for i=1:modelnr
                fullpeak=fullpeak+mod_ga(gx,parar(:,i));
            end
            
            [ma,id]=max(fullpeak);
            left=find(fullpeak>ma/300,1,'first');
            right=find(fullpeak>ma/300,1,'last');
            if info
                disp('Split Peak Right        ');
                disp(['    Left bound:        ' num2str(gx(left))]);
                disp(['    Peak max location: ' num2str(gx(id))]);
                disp(['    Right bound:       ' num2str(gx(right))]);
                disp(['    Peak height:       ' num2str(ma)]);
                disp(['    Peak rel. height:  ' num2str(ma/ma_global)]);
                disp(['    Model number:      ' num2str(modelnr)]);
            end
            p.peakX_l2r=gx(left:right);
            p.fullpeak_l2r=fullpeak(left:right);
            p.paramat=parar;
            p.lbound=gx(left);
            p.rbound=gx(right);
            p.maxlocation=gx(id);
            p.height=ma;
            p.relheight=ma/ma_global;
            p.modelnumber=modelnr;
            
            pr=peak(p.peakX_l2r,p.fullpeak_l2r,modelnr,parar,p);
            
        end
        
        %% UPDATE PEAK INFO AFTER E.G. CORRECTION HAS BEEN CHANGED
        % !! USE IF CORR IS INSIDE OF LOCAL GRID
        function obj=updatePropAfterCorr(obj,gx,ma_global)
            % Inputs:
            % gx - full time grid
            % ma_global - signal maximum on full grid
            fullpeak=obj.getmod(gx);
            [ma,id]=max(fullpeak);
            left=find(fullpeak>ma/300,1,'first');
            right=find(fullpeak>ma/300,1,'last');
            obj.peakstruct.peakX_l2r=gx(left:right);
            obj.peakstruct.fullpeak_l2r=fullpeak(left:right);
            obj.peakstruct.lbound=gx(left);
            obj.peakstruct.rbound=gx(right);
            obj.peakstruct.maxlocation=gx(id);
            obj.peakstruct.height=ma;
            obj.peakstruct.relheight=ma/ma_global;
            
            cor=interp1(obj.timegrid,obj.correction,gx,'linear',0);
            obj.correction=cor(left:right);
            
            obj.timegrid=gx(left:right);
            obj.signal=fullpeak(left:right);
        end
        
        %% ADD A CORRECTION 
        % !! USE IF CORR IS OUTSIDE OF LOCAL GRID
        function obj=addCorrection(obj,gx,ma_global,corr)
            % Inputs:
            % gx - full time grid
            % ma_global - signal maximum on full grid
            % corr - correction on full grid to process
            
            fullpeak=obj.getmod(gx); % get model + current correction
            fullpeak=fullpeak+corr; % add new correction
            [ma,id]=max(fullpeak);
            left=find(fullpeak>ma/300,1,'first');
            right=find(fullpeak>ma/300,1,'last');
            obj.peakstruct.peakX_l2r=gx(left:right);
            obj.peakstruct.fullpeak_l2r=fullpeak(left:right);
            obj.peakstruct.lbound=gx(left);
            obj.peakstruct.rbound=gx(right);
            obj.peakstruct.maxlocation=gx(id);
            obj.peakstruct.height=ma;
            obj.peakstruct.relheight=ma/ma_global;
            
            fullCorrOld=interp1(obj.timegrid,obj.correction,gx,'linear',0);
            obj.correction=corr(left:right)+fullCorrOld(left:right);
            
            obj.timegrid=gx(left:right);
            obj.signal=fullpeak(left:right);
        end
        
        %% PLOT
        function p=plot(obj,col,lw,ax)
            if nargin<4
                ax=gca;
            end
            p=[];
            if nargin==1
                p=plot(ax,obj.timegrid,obj.getmod(obj.timegrid),'k');
            elseif nargin==2
                p=plot(ax,obj.timegrid,obj.getmod(obj.timegrid),'color',col);
            elseif nargin==3 || nargin==4
                p=plot(ax,obj.timegrid,obj.getmod(obj.timegrid),'color',col,'Linewidth',lw);
            end
        end
        
        %% FILL PLOT
        function p=fill(obj,col,ax)
            if nargin<3
                ax=gca;
            end
            p=fill(ax,obj.timegrid,obj.getmod(obj.timegrid),col,'facealpha',0.7);
            
        end
        
        %% DO REFIT OF THE PEAK
        function pnew=refit(obj,fullgx,fullgc,fullmod)
            % Output: 
            % a new peak with the same sub peak number 
            % the previous correction is ignored 
            % -> new peak has zero correction
            
            % debug plots?
            isdebug=0;
            
            %get bounds
            lb=obj.peakstruct.lbound;
            rb=obj.peakstruct.rbound;
            
            %find indices in full profile
            lbx=find(fullgx<=lb,1,'last');
            rbx=find(fullgx>=rb,1,'first');
            
            %get gc and model window
            gxt=fullgx(lbx:rbx);
            gct=fullgc(lbx:rbx);
            modt=fullmod(lbx:rbx);
            
            %get residuals without the current peak
            rest=gct-modt+obj.getmod(gxt);
            
            %Scaling to max 1 / get center idx
            [scal1,idce]=max(rest);
            rest_scal=rest/scal1;
            center=gxt(idce);
            
            %center to max position
            gxt_shift=gxt-center;
            
%             % approximated width
%             widthL=0.15*(rb-lb);
%             widthR=0.15*(rb-lb);
%             
%             hiscale=0.9;
            
            % get width by approximate when the peak reaches 3/5 height
            idxL=find(((gxt_shift<=0)+(   rest_scal<=3/5   ))==2,1,'last');
            widthL=-gxt_shift(idxL);
            idxR=find(((gxt_shift>=0)+(   rest_scal<=3/5   ))==2,1,'first');
            widthR=gxt_shift(idxR);
            
            if isempty(widthL) || isempty(widthR)
                widthid=1;
                disp(['<> getApproxPeak: widthL/R empty'])
            else
                if widthL<widthR
                    xid=idxL;
                    width=widthL;
                else
                    xid=idxR;
                    width=widthR;
                end
                widthid=abs(idce-xid);
            end
            
            % tight peak range (constraint model=signal)
            le_idx=max([1 idce-floor(0.5*widthid)]);
            ri_idx=min([length(rest_scal) idce+ceil(0.5*widthid)]);
            
            %active windows
            actx=true(size(gxt_shift));
%             actx(le_idx:ri_idx)=true;
            actx_w=true(size(gxt_shift));
            
            
% %             if obj.isdebug
% %                 ft=figure;
% %                 subplot(2,2,1)
% %                 plot(gxt,gct,'k',gxt, modt,'b','linewidth',2);
% %                 legend('signal','model');
% %                 yy=ylim;
% %                 
% %                 subplot(2,2,2)
% %                 plot(gxt_shift,rest_scal,'k','linewidth',2);
% %                 hold on
% %                 plot(gxt_shift([1 end]),obj.ga.rescl.noiselvl/scal1*[1 1],'r--','linewidth',2)
% %                 plot(0,hiscale,'g*','linewidth',2);
% %                 plot([-widthL widthR],0.5*[1 1],'g-*','linewidth',2);
% %                 plot(gxt_shift(actx),zeros(size(gxt_shift(actx))),'m-','linewidth',2)
% %                 legend('residual','noise level','max','width','data fit win');
% %                 title('shifted/scaled')
% %             end
            
            
            % set initial values and lower/upper bounds
            hiscale=0.95;
            
            p0=obj.parameters;
            p0(1,:)=p0(1,:)-center; %shift centers
            p0(4,:)=hiscale*(p0(4,:)./scal1); %scale height
            
            
           
            pL=zeros(size(p0));
            pR=inf(size(p0));
            pL(1)=min(gxt_shift);
            pR(1)=max(gxt_shift);
            pL(2)=0.05*p0(2);
            pL(3)=0.05*p0(3);
            pR(2)=2.5*p0(2);
            pR(3)=10*p0(3);
            pR(5)=1;
            pR(6)=1;
            pL=pL(:);
            pR=pR(:);
            p0=p0(:);
            
            
            
            %fit
            we=[1 0];
            fun2=@(p) tar_opt(p,gxt_shift,rest_scal,size(obj.parameters,2),6,we,[],[],actx,actx_w,0);
            options=optimoptions(@lsqnonlin,...
                'Display','iter',...
                'MaxIter',100,...
                'MaxFunEvals',10000,...
                'TolFun',1e-10,...
                'TolX',1e-10);
            pOpt=lsqnonlin(fun2,p0,pL,pR,options);
            
            pOpt=reshape(pOpt,6,size(obj.parameters,2));
            
            %model eval / rescale / reshift
            mod_gxt_shift=zeros(length(gxt_shift),1);
            for i=1:size(pOpt,2)
                mod_gxt_shift=mod_gxt_shift+mod_ga(gxt_shift,pOpt(:,i));
            end
            pOpt(4,:)=scal1*pOpt(4,:);
            pOpt(1,:)=pOpt(1,:)+center;
            
            
            
            fullpeak=zeros(length(fullgx),1);
            for i=1:size(pOpt,2)
                fullpeak=fullpeak+mod_ga(fullgx,pOpt(:,i));
            end
            
            if isdebug
                figure();
                subplot(2,2,3)
                plot(gxt_shift,rest_scal,'k','linewidth',2);
                hold on
                plot(gxt_shift,mod_gxt_shift,'r--','linewidth',2);
                legend('residual','model');
                title('shifted/scaled')
                xx=xlim;
                
                subplot(2,2,4)
                plot(fullgx,fullgc,'k','linewidth',2);
                hold on
                plot(fullgx,fullpeak,'r--','linewidth',2);
                xlim(xx+[-(xx(2)-xx(1)) (xx(2)-xx(1))]+center)
            end
            
            %build peak
            ma_global=max(fullgc);
            [ma,id]=max(fullpeak);
            left=find(fullpeak>ma/300,1,'first');
            right=find(fullpeak>ma/300,1,'last');
            if isempty(left)
                left=1;
            end
            if isempty(right)
                right=length(fullpeak);
            end
            inte=trapz(fullgx,fullpeak);
            pp.peakX_l2r=fullgx(left:right);
            pp.fullpeak_l2r=fullpeak(left:right);
            pp.paramat=pOpt;
            pp.lbound=fullgx(left);
            pp.rbound=fullgx(right);
            pp.maxlocation=fullgx(id);
            pp.height=ma;
            pp.relheight=ma/ma_global;
            pp.modelnumber=1;
            
            pnew=peak(pp.peakX_l2r,pp.fullpeak_l2r,pp.modelnumber,pOpt,pp);
        end
        
        
    end
end














