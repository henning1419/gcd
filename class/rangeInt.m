classdef rangeInt < profileC
    %RANGEINT (see peak class for details)
    
    properties
        %signal
        %timegrid
        
        lbound
        rbound
        height
        maxlocation
        intinfo=[];
    end
    
    methods
        %% Constructor
        function obj = rangeInt()
        end
        
        %% initialize
        function obj = init(obj,bounds,gx,gc,area_global)
            if isempty(bounds)
                return
            end
            
            validateattributes(bounds,{'double'}, {'numel', 2,'increasing'})
            
            obj.lbound=bounds(1);
            obj.rbound=bounds(2);
            
            idxRangeBool=false(size(gx));
            idxRangeBool(gx>=obj.lbound & gx<=obj.rbound)=true;
            
            % get segment of signal and grid
            obj.timegrid=gx(idxRangeBool);
            obj.signal=gc(idxRangeBool);
            
            [obj.height,idx]=max(obj.signal);
            obj.maxlocation=obj.timegrid(idx);
            
            obj.intinfo.sigAbs=trapz(obj.timegrid,obj.signal);
            obj.intinfo.sigRel=obj.intinfo.sigAbs/area_global;
        end
        
        
    end
end

