classdef (Abstract) hasHistory < handle
    %HASHISTORY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(GetAccess = public, SetAccess = protected)
        isHistory=1;
        history={};
    end
    
    
    methods
        
        %% ENABLE HISTORY
        function enableHistory(obj)
            obj.isHistory=1;
        end
        
        %% DISABLE HISTORY
        function disableHistory(obj)
            obj.isHistory=0;
        end
        
        %% ADD HISTORY ENTRY
        function addHist(obj,comment)
            if ~obj.isHistory
                return;
            end                        
            st=comment;
            st.timestamp = datestr(now,'yyyy-mm-dd HH:MM:SS');
            
            obj.history{end+1}=st;
        end
        
        %% DELETE LAST HISTORY ENTRY
        function delLastHist(obj)
            if ~isempty(obj.history)
                obj.history{end}=[];
            end
        end
        
        %% PRINT HISTORY
        function printHistory(obj,fileID)
            if nargin==1
                isFileOutput=0;
            else
                isFileOutput=1;
            end
            
            if isFileOutput
                fprintf(fileID,'History info:\n\n');
                
                fprintf(fileID,'%10s %10s %10s %12s %10s %10s %10s\n','Timestamp','Cluster idx', 'Cluster active','Job success','Job idx');% ,'calcWindow', 'Deadzone');
                for i=1:length(obj.history)
                    
                    if obj.history{i}.success == -1 % job not started
                        curr_suc = 'not started';
                    elseif obj.history{i}.success == 0 % job failed
                        curr_suc = 'failed';
                    else % job succeeded
                        curr_suc = 'succeeded';
                    end
                   
                    fprintf(fileID,'%10f %10f %10f %12s %10f\n',obj.history{i}.timestamp, obj.history{i}.cluster, obj.history{i}.clusterAct, curr_suc, obj.history{i}.jobidx);%,obj.history{i}.calcWindow,obj.history{i}.deadzone);
                    
                    
%                     hist_comment.jobidx = jobidx;
%                     hist_comment.deadzone = obj.clusterP(ClusterIdx).deadzone;
%                     hist_comment.success = -1; % job not started
%                     hist_comment.cluster = ClusterIdx;
%                     hist_comment.calcWindow = clusterCurCalcZone{ClusterIdx};
%                     hist_comment.clusterAct = clusterAct(ClusterIdx);
%                     hist_comment.pOpt = resJob.pOpt;
%                     hist_comment.p = p;
                end
                
                fprintf(fileID,'\n');
                
            else
                fprintf('History info:\n\n');
                
                fprintf('%10s %10s %10s %12s %10s %10s %10s\n','Timestamp','Cluster idx', 'Cluster active','Job idx','Job success');% ,'calcWindow', 'Deadzone');
                for i=1:length(obj.history)
                    
                    if obj.history{i}.success == -1 % job not started
                        curr_suc = 'not started';
                    elseif obj.history{i}.success == 0 % job failed
                        curr_suc = 'failed';
                    else % job succeeded
                        curr_suc = 'succeeded';
                    end
                   
                    fprintf('%10f %10f %10f %12s %10f \n',obj.history{i}.timestamp, obj.history{i}.cluster, obj.history{i}.clusterAct, curr_suc, obj.history{i}.jobidx);%, obj.history{i}.calcWindow,obj.history{i}.deadzone);
                end
            end
            fprintf('\n')
        end
    end
end

