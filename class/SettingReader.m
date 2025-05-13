classdef SettingReader < handle
    
    properties
        % DEFAULT PROPERTY VALUES
        
        info    = 1;      % Info-Details (0=aus, 1=normal, 2=alles)
        ver
        d       = 3;      % SG-Filter - Polynomgrad
        M       = 10;     % SG-Filter - halbe Fensterbreite
        
        bl_n_all
        bl_relmaxpeakW      = 0.43;
        bl_PP_relmaxpeakW   = 1/2;
        bl_startend         = 120;
        bl_PPruns           = 4;
        
        cl_n                = 300;      % Cluster - #Gitterpunkte (500 für gesamtes GC-Profil)
        cl_first            = 100;      % Cluster - Noiselevel mittels kleinsten ... Pkt.
        cl_fac              = 1.3;      % Cluster - Skalierung für Noiselevel 1 - 2
        NoiseWin            = 501;      % Cluster - MovingWindow width in Noiselevel approx (has to be odd, >301)
        
        maxpeak_nr                      % Fit - Maximale Anzahl an Peaks je Subproblem
        fev_tol             = 0.3;      % Fit - Quadratsum.-toleranz der Optimierung, Peak verwerfen?
        errTol              = 0.001;    % Fit - Abbruch bei errTol Restintegral 0.001 fein, 0.01 grob
        retry               = 50;
        epsi                = 0.001;
        minPdist            = 0.01;
        delta               = -1e-2;
        peakoverlap         = 0.01;
        overlap_runs        = 10;
        relPeakHi           = 1/6;
        relValley           = 1/2;
        relValleyBoth       = 4/5;
        
        ispar               = 1;
        numWorkers
        pooloverload        = 0;
        isDetailedPlot      = 0;
        maxCalcPerClus      = 1;
        bw
        isREF
        dataR
        tempf
        
    end
    
    methods
        %% Constructor
        function obj = SettingReader() 
        end
        
        %% Program-specific parameters are set
        function setParameter(obj, paraini)
            obj.ver = paraini.ver;
            obj.bl_n_all = paraini.bl_n_all;
            obj.maxpeak_nr = paraini.maxpeak_nr;
            obj.bw = paraini.bw;
            obj.isREF = paraini.isREF;
            if paraini.isREF
                obj.dataR = paraini.dataR;
            end
            obj.tempf = paraini.tempf;
            obj.numWorkers = max(1,feature('numcores')-1);
        end

        %% Validation
        function validation(obj,path_para)
            % check if the file type is correct
            [~,~,ext] = fileparts(path_para);
            zulEndung=['.txt','.dat','.csv','.xls','.xlsb','.xlsm', '.xlsx', '.xltm', '.xltx', '.ods', 'xml'];
            if ~any(ismember(ext, zulEndung))
                error('The parameter file has an invalid file type and is not considered.')
            end            
            
            T = readtable(path_para);
            p = properties(obj);
            
            % check if file is invalid due to the dimension of T
            if size(T,1)>length(p)
                warning('Too many parameters have been specified, some will be ignored.')
            end
            if size(T,2)~=2
                error('Wrong Parameter specification (too many or too few values assigned). Parameter file is not considered.')
            end
            
            % check if all names in the parameter file are properties
            T.Properties.VariableNames{1}='Var1';
            for i=1:length(T.Var1)
                if isempty(find(strcmp(T.Var1{i},p),1))
                    warning([T.Var1{i}, ' is not a property and will be ignored.'])
                end
            end
        end
        
        %% Read parameter file
        function para = read(obj,path_para) 
            % input: path to settings txt-file, see ./bin/parameter_settings.txt for an example
            % output: para struct
            
            if ~isempty(path_para)
                %read from file
                
                % validate parameter file
                obj.validation(path_para)
                
                
                T = readtable(path_para);
                p = properties(obj);
                T.Properties.VariableNames{1}='Var1';
                T.Properties.VariableNames{2}='Var2';
                
                % parameter values are checked and overwrite class default
                % properties
                fprintf('%20s %14s %14s\n','Parameter name','Value','Source')
                for i=1:length(p)
                    idx = find(strcmp(p{i},T.Var1),1,'first');
                    if  any(idx)
                        if ~isa(T.Var2(idx(1)),'numeric') || isnan(T.Var2(idx(1)))
                            error(['Value of Parameter ', num2str(p{i}),' is invalid. Parameter file is not considered.'])
                        end
                        obj.(p{i}) = T.Var2(idx(1));
                        src='external';
                    else
                        src='default';
                    end
                    if isnumeric(obj.(p{i}))
                        fprintf('%20s %14d %14s\n',p{i},obj.(p{i}),src)
                    else % string
                        fprintf('%20s %14s %14s\n',p{i},obj.(p{i}),src)
                    end
                end
            else
                %only defaults
                fprintf('%20s %14s %14s\n','Parameter name','Value','Source')
                src='default';
                p = properties(obj);
                for i=1:length(p)
                    if isnumeric(obj.(p{i}))
                        if isscalar(obj.(p{i}))
                            fprintf('%20s %14d %14s\n',p{i},obj.(p{i}),src)
                        else
                            fprintf('%20s %14s %14s\n',p{i},'vector',src)
                        end
                    else % string
                        fprintf('%20s %14s %14s\n',p{i},obj.(p{i}),src)
                    end
                end
            end
            
            % write settings to para struct
            for i=1:length(p)
                eval(['para.' p{i} '=obj.' p{i} ';' ]);
            end
        end
            
    end
end
    
