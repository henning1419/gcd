classdef Val2ValComp
    
    properties (Access=public)
        val
        inf
    end
    
    methods (Access=public)
        function obj = Val2ValComp(val,inf)
            %constructor
            obj.val=val;
            obj.inf=inf;
            
        end
        
        function obj = plot(obj)
            for i=1:length(obj.val)
                obj.val(i).plot;
                set(gcf,'Name',obj.inf{i});
            end
        end
        
        function obj = time_avarage(obj)
            disp(' ')
            disp(['======> Time Average'])
            for i=1:length(obj.val)
                disp(['### ' obj.inf{i} ' ###'])
                obj.val(i).time_average;
            end
        end
        
        function obj = time_total(obj)
            disp(' ')
            disp(['======> Time Total'])
            for i=1:length(obj.val)
                disp(['### ' obj.inf{i} ' ###'])
                obj.val(i).time_total;
            end
        end
        
        function obj = range_plot(obj,id,range)
            if nargin==2
                range=[];
            end
            
            ax_arr=[];
            figure('Position',[100 100 1000 800]);
            if ~isempty(range)
                set(gcf,'Name',['ID: ' num2str(id) ' - Range: ' num2str(range(1)) ' - ' num2str(range(2))]);
            else
                set(gcf,'Name',['ID: ' num2str(id)]);
            end
            for i=1:length(obj.val)
                ga=obj.val(i).loadGA(id);
                ax=subplot(length(obj.val),1,i);
                ax_arr(end+1)=ax;
                ga.profPF.plot('k',2,ax);
                hold on
                p=ga.modePF.plot('r',2,ax);
                set(p,'Linestyle','--')
                title([obj.inf{i} ' (n=' num2str(length(ga.profPF.timegrid)) ')'],'Interpreter','none')
                xlabel('Retentionszeit')
                ylabel('Intensitaet')
                
                %view area
                if ~isempty(range)
                    xlim(range)
                    gx=ga.profPF.timegrid;
                    idL=find(gx<range(1),1,'last');
                    idR=find(gx>range(2),1,'first');
                    ma=max(ga.profPF.signal(idL:idR));
                    ylim([min(ga.profPF.signal(idL:idR)) 1.05*ma]);
                end
            end
            linkaxes(ax_arr)
            
        end
    end
end
