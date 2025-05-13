classdef profileC < handle
    %%
    properties
        % timegrid - double vector - retention time grid 
        timegrid 
        
        % signal - double vector
        signal
        
    end
    
    %% methods
    methods
        %% constructor
        function obj = profileC(timegrid,signal)
            if nargin==2
                obj.timegrid = timegrid;
                obj.signal =signal;
            elseif nargin==0
                
            end
        end
        
%         %% destructor
%         % Destructor
%         function delete(obj)
%             p = properties(obj);
%             for i=1:length(p)
%                 try 
%                     eval(['delete(obj.' p{i} ')'])
%                 catch
%                     eval(['obj.' p{i} '=[];'])
%                 end
%             end
%             clear obj
%         end
        
        %% baseline (deprecated, see GCanalysis.getBaseline)
        function resbl = calcBaselineP(obj,para,ng)
            d=para.d;
            M=para.M;
            if nargin==2
                bl_n_all=para.bl_n_all;
            elseif nargin==3
                bl_n_all=ng;
            end
            
            resbl=corr_bl(obj.timegrid,...
                obj.signal,...
                bl_n_all,d,2*M+1);
        end
        
        %% baseline_testing (deprecated, see GCanalysis.getBaseline)
        function resbl = calcBaselineP2(obj,para,ng)
            d=para.d;
            M=para.M;
            if nargin==2
                bl_n_all=para.bl_n_all;
            elseif nargin==3
                bl_n_all=ng;
            end
            
            resbl=corr_bl(obj.timegrid,...
                obj.signal,...
                bl_n_all,d,2*M+1);
        end
        
        %% clustering (used in GCanalysis.calcCluster)
        function rescl = calcClusterP(obj,para)
            
            rescl=corr_cluster(obj.timegrid,...
                obj.signal,...
                para);
        end
        
        
        %% plot
        function p=plot(obj,col,lw,ax)
            if nargin<4
                ax=gca;
            end
            p=[];
            if nargin==1
                p=plot(ax,obj.timegrid,obj.signal,'k');
            elseif nargin==2
                p=plot(ax,obj.timegrid,obj.signal,'color',col);
            elseif nargin==3 || nargin==4
                p=plot(ax,obj.timegrid,obj.signal,'color',col,'Linewidth',lw);
            end
        end  
        
    end
 
end

