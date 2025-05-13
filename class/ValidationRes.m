classdef ValidationRes
    %VALIDATIONRES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access=public)
        mat_list     % file list of matlab result files
        labels
        time_vals
        accuracy_vals
        ds_length
        ds_nr_peaks
        len
        ref_arr
    end
    
    methods (Access=public)
        function obj = ValidationRes(mainfolder,isRef)
            %constructor
            
            if nargin==1
                isRef=0;
            end
            
            addpath('./bin')
            cpath=cd;
            cd(mainfolder);
            templist=dir('**/*.mat');
            cd(cpath)
            obj.mat_list={};
            if isRef
                obj.ref_arr={};
            end
            obj.labels={};
            obj.time_vals=[];
            obj.accuracy_vals=[];
            obj.ds_length=[];
            obj.ds_nr_peaks=[];
            h=waitbar(0,['loading result data 0 of ' num2str(length(templist))]);
            obj.len=length(templist);
            for i=1:obj.len
                tempf=[templist(i).folder filesep templist(i).name];
                disp(['Loading from: ' tempf])
                obj.mat_list{end+1}=tempf;
                obj.labels{end+1}=templist(i).name;
                load(tempf,'ga','para');
                obj.time_vals(end+1)=ga.progress(end);
                if isRef
                    obj.ref_arr{end+1}=para.dataR;
                end
                obj.accuracy_vals(end+1)=ga.coveredPeakAreaNoN;
                obj.ds_length(end+1)=length(ga.baseP.timegrid);
                obj.ds_nr_peaks(end+1)=length(ga.peaks);
                clear ga
                clear para
                
                if ~ishandle(h)
                    h=waitbar(i/length(templist),['loading result data ' num2str(i) ' of ' num2str(length(templist))]);
                else
                    waitbar(i/length(templist),h,['loading result data ' num2str(i) ' of ' num2str(length(templist))]);
                end
            end
            close(h)
        end
        
        function r = time_total(obj)
            r=sum(obj.time_vals);
            tl='s';
            if r>60
                r=r/60;
                tl='min';
            end
            if r>60
                r=r/60;
                tl='h';
            end
            if nargout==0
                disp(['Total computation time for ' num2str(length(obj.time_vals)) ' profiles: ' num2str(r) ' ' tl]);
                disp(' ')
            end
        end
        
        function r = time_average(obj)
            r=sum(obj.time_vals)/length(obj.time_vals);
            mi=min(obj.time_vals);
            ma=max(obj.time_vals);
            tl='s';
            if r>60
                r=r/60;
                mi=mi/60;
                ma=ma/60;
                tl='min';
            end
            if r>60
                r=r/60;
                mi=mi/60;
                ma=ma/60;
                tl='h';
            end
            if nargout==0
                disp(['Average computation time: ' num2str(r) ' ' tl]);
                disp(['Minimum computation time: ' num2str(mi) ' ' tl]);
                disp(['Maximum computation time: ' num2str(ma) ' ' tl]);
                disp(' ')
            end
        end
        
        
        
        function plot(obj,info)
            if nargin==1
               info=''; 
            end
            figure()
            subplot(2,1,1)
            histogram(obj.time_vals./60,20);
            xlabel('computation time [min]')
            ylabel('number of data sets')
            title(info)
            
            subplot(2,1,2)
            histogram(obj.accuracy_vals,20);
            xlabel('accuracy (noise-corrected)')
            ylabel('number of data sets')
        end
        
        function plotTimePerLength(obj)
            figure()
            plot(obj.ds_length,obj.time_vals./60,'*');
            
            xlabel('gc data points')
            ylabel('comp time [min]')
        end
        
        function plotTimePerPeak(obj,info)
            if nargin==1
               info=''; 
            end
            figure()
            plot(obj.ds_nr_peaks,obj.time_vals./60,'*');
            title(info)
            xlabel('number of extracted peaks')
            ylabel('comp time [min]')
        end
        
        function plot_single(obj)
            figure()
            subplot(2,1,1)
            for id=1:obj.len
                h=plot(obj.time_vals(id)./60,0,'ro');
                set(h,'ButtonDownFcn',{@fcn_label,obj,id})
                hold on
            end
            xlabel('computation time [min]')
            
            subplot(2,1,2)
            for id=1:obj.len
                h=plot(obj.accuracy_vals(id),0,'ro');
                set(h,'ButtonDownFcn',{@fcn_label,obj,id})
                hold on
            end
            xlabel('accuracy (noise-corrected)')
        end
        
        function ga=loadGA(obj,id)
            cf=obj.mat_list(id);
            load(cf{1},'ga')
        end
        
        function plotGA(obj,id)
            cf=obj.mat_list(id);
            load(cf{1},'ga')
            plot(ga)
        end
        
        function plotGApeaks(obj,id)
            cf=obj.mat_list(id);
            load(cf{1},'ga')
            plotPeaks(ga)
        end
        
        function plotGABL(obj,id)
            cf=obj.mat_list(id);
            load(cf{1},'ga')
            plotBL(ga)
        end
        
        function [ax1,ax2]=plotGARef(obj,id)
            cf=obj.mat_list(id);
            load(cf{1},'ga')
            [ax1,ax2]=plotRefComp(ga);
        end
        
        function printfRange(obj,range,fol)
            h=waitbar(0,['saving picture 0 of ' num2str(obj.len)]);
            for i=1:obj.len
                if ~ishandle(h)
                    h=waitbar(i/obj.len,['saving picture ' num2str(i) ' of ' num2str(obj.len)]);
                else
                    waitbar(i/obj.len,h,['saving picture ' num2str(i) ' of ' num2str(obj.len)]);
                end
                
                cf=obj.mat_list(i);
                load(cf{1},'ga')
                f=ga.plotSingle(0);
                
                %view area
                xlim(range)
                gx=ga.profPF.timegrid;
                idL=find(gx<range(1),1,'last');
                idR=find(gx>range(2),1,'first');
                ma=max(ga.profPF.signal(idL:idR));
                ylim([min(ga.profPF.signal(idL:idR)) 1.05*ma]);
                
                %Save pic
                print('-dpng','-loose',[fol filesep 'range_' num2str(range(1)) ...
                    '_' num2str(range(2)) '_' num2str(i) '.png']);
                delete(f)
                clear ga
            end
        end
        
        function printfMax(obj,fol,relWinWidth)
            if nargin~=3
                relWinWidth=0.1;
            end
            
            h=waitbar(0,['saving picture 0 of ' num2str(obj.len)]);
            for i=1:obj.len
                if ~ishandle(h)
                    h=waitbar(i/obj.len,['saving picture ' num2str(i) ' of ' num2str(obj.len)]);
                else
                    waitbar(i/obj.len,h,['saving picture ' num2str(i) ' of ' num2str(obj.len)]);
                end
                
                cf=obj.mat_list(i);
                load(cf{1},'ga')
                f=ga.plotSingle(0);
                
                %view area
                gx=ga.profPF.timegrid;
                gc=ga.profPF.signal;
                [ma,maidx]=max(gc);
                lengc=length(gc);
                lenwin=round(lengc*relWinWidth/2);
                
                idL=max(1,maidx-lenwin);
                idR=min(lengc,maidx+lenwin);
                range=[idL idR];
                
                xlim(gx(range))
                ylim([min(gc(idL:idR)) 1.05*ma]);
                
                %Save pic
                print('-dpng','-loose',[fol filesep 'max_' num2str(i) '.png']);
                delete(f)
                clear ga
            end
            close(h)
        end
        
        function printMatList(obj)
            for i=1:obj.len
                disp([num2str(i) ': ' obj.mat_list{i}]);
            end
        end
        
        function printRefAnz(obj,relHi)
            
            ref_nr_peaks=zeros(obj.len,1);
            for i=1:obj.len
                cref=obj.ref_arr{i};
            
                ch=zeros(1,size(cref,1));
                
                for j=1:size(cref,1)
                    ch(j)=cref(j,4);
                end 
                chmax=max(ch);
                
                for j=1:length(ch)
                    if ch(j)/chmax>=relHi
                        ref_nr_peaks(i)=ref_nr_peaks(i)+1;
                    end
                end 
            
            end
            
            
            
%             for i=1:obj.len
%                 
%                 ref_nr_peaks(i)=size(cref,1);
%             end 
            
            mod_nr_peaks=zeros(length(obj.ds_length),1);
            h=waitbar(0,['processing 0 of ' num2str(obj.len)]);
            for i=1:obj.len
                if ~ishandle(h)
                    h=waitbar(i/obj.len,['processing ' num2str(i) ' of ' num2str(obj.len)]);
                else
                    waitbar(i/obj.len,h,['processing ' num2str(i) ' of ' num2str(obj.len)]);
                end
                
                cf=obj.mat_list(i);
                load(cf{1},'ga')
                
                ch=zeros(1,length(ga.peaks));
                
                for j=1:length(ga.peaks)
                    ch(j)=ga.peaks(j).peakstruct.height;
                end 
                chmax=max(ch);
                
                for j=1:length(ch)
                    if ch(j)/chmax>=relHi
                        mod_nr_peaks(i)=mod_nr_peaks(i)+1;
                    end
                end 

                clear ga
            end
            close(h)

            figure
            plot(1:obj.len, mod_nr_peaks(:)-ref_nr_peaks(:),'rx','Linewidth',2,'Markersize',12);
            xlim([1 obj.len])
            
            idx=find(mod_nr_peaks(:)-ref_nr_peaks(:)<0);
            disp(['data sets with more ref peaks than gca peaks for relHi ' num2str(relHi)])
            disp(idx)
            
        end
        
        
    end
end

