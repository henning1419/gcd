function  main(tempf,mode,varargin)
%% Application cases
% win utf8 csv
% [win utf16 csv via knime (utf8 conversion)]
% win ch
% unix utf16 csv
% [unix utf8 csv (modify load_gc_singleFile)]
% unix ch

% recommendated format: ch

%% INPUT
% all inputs are strings

% tempf        - full filepath with gc data (col 1: grid, col 2: signal) or
%                of result mat-file (depends on mode)
% mode         - {calc, analyze, plot}
%                calc    - data load, preprocess, fitting, postprocess for gc data
%                          in tempf
%                analyze - start GUI for a result mat-file in tempf
%                new     - start GUI without prior analysis. loading is
%                          first step
%                plot    - make a simple evaluation plot for a result mat-file in tempf
% exclude      - format: name-value-pair ->  'exclude','[3,5;20,60;70,72]' 
% (optional)     (matlab matrix-like, 2 columns)
%                each row defines a range of retention times
%                retriction: row-wise strictly ascending values, e.g. 3<5, 20<60, 70<72
%                remarks: overlaping ranges are concatenated, validation of
%                restrictions will throw an error
% intRange     - format: name-value-pair ->  'intRange','[3,5;20,60;70,72]' 
% (optional)     (matlab matrix-like, 2 columns)
%                each row defines a range of retention times
%                retriction: row-wise strictly ascending values, e.g. 3<5, 20<60, 70<72
%                remarks: validation of restrictions will throw an error
% settingsPath - format: name-value-pair ->  'settingsPath','/path2settings/file.txt'
% (optional)     remarks: an example file can be found in the subfolder /bin

%% CALL EXAMPLES
% main([],'calc')
% main(path,'calc')
% main([],'new')
% main(path,'analyze')
% main(path,'calc','exclude','[20,50;90,130]')
% main(path,'calc','exclude','[20,50;90,130]','intRange','[70,80;140,150]')
% main(path,'calc','settingsPath','/path2settings/file.txt')


%% REMARKS
% exlude ranges have to be considered in all error evaluations 
% and most fitting processes.
% peak.refit is only used in GUI analyzer, the exclude ranges are not
% considered here.


%% TODO
% signal shift by noise (pre)
% modify utf8 conversion script from "apply to folder" to "apply to single file"
% split correction
% apply "adjust tailing peaks" to consecutive tailing peaks (currently only pairs of tailing peaks)

ver='3.0.0';


%% START CODE

disp(['GC Analysis - Version ' ver])

%% PROCESS INPUTS
def_exRange='[]';
def_inRange='[]';
def_inSettingsPath=[];
expectedModes = {'calc','plot','analyze','new'};

validationFcnRange = @(x) validateattributes(x,{'char'},{'nonempty'});

pa=inputParser;
pa.addRequired('tempf');
pa.addRequired('mode',@(x) any(validatestring(x,expectedModes)));
pa.addParameter('exclude',def_exRange,validationFcnRange);
pa.addParameter('intRange',def_inRange,validationFcnRange);
pa.addParameter('settingsPath',def_inSettingsPath);

pa.parse(tempf,mode,varargin{:});

exRange=pa.Results.exclude;
intRange=pa.Results.intRange;
settingsPath=pa.Results.settingsPath;


       

%% MODE SELECTION
switch mode
    case 'calc'
        
        isplot=1;
        
        
        addpath('./class')
        addpath('./bin')
        
        % tempf is used as path+file
        [gc,gx,ident,flag]=load_gc_csv_singleFile(tempf);
        if flag==0
            error('loading failed')
        end

        % adaptive thinning
        % points per second
        pps=length(gx)/(gx(end)-gx(1));
        tar_pps=400;
        if tar_pps<pps
            nr=ceil((gx(end)-gx(1))*tar_pps);
            gxn=linspace(gx(1),gx(end),nr);
            gcn=interp1(gx,gc,gxn);

            gc=gcn;
            gx=gxn;
        else
            warning('Target pps is larger than raw pps. Low resolution device/measurement ?');
        end
        
        ident=['lastCalc'];
        isREF=0;
              
        savefol=[cd filesep ident filesep];
        
        gx=gx(:);
        gc=gc(:);
        gc=gc-min(gc);
        
        %% Setting parameters
        
        % initialization with default parameters
        % para contains the parameters as public class properties
        % (alternatively para can be defined as struct)
        sr = SettingReader();
        
        % define and add programm specific parameters
        lenref=ceil(gx(end)-gx(1));
        paraini.ver = ver;
        paraini.bl_n_all = max(ceil(1/10*lenref),10);
        paraini.maxpeak_nr = max(ceil(12*lenref),400);
        paraini.bw = find((gx-gx(1))>sr.minPdist,1,'first');
        paraini.isREF = isREF;
        if isREF
            paraini.dataR = dataR;
        end
        paraini.tempf = savefol;
        
        sr.setParameter(paraini);
        
        % due to given parameter file the parameters will be updated
        para=sr.read(settingsPath);
        
%         para.numWorkers=1;
        %% START
        % Info: Wenn hier Fehler kommen kann es daran liegen, dass das
        % adaptive thinning die Daten "kaputt" macht + auf nur einem Worker
        
        
        mkdir(savefol)
                
        if exist('ga','var')
            clear('ga')
        end
        
        ga=GCanalysis(gx,gc,para);
        ga.startParPool;
        
        % Add file path+name to ga object
        ga.para.fName = tempf;
        
        ga.loadExcludeRanges(exRange);
        ga.loadIntRanges(intRange);
        
        ga.calcBaseline;
        ga.calcCluster;
        %update cluster deadzones according to exclude ranges
        ga.updateClusterDeadzones;
        ga.dispProgress;
        
        ga.enableHistory;
        
        ga.calcFitPar;
        ga.dispProgress;
        
        ga.generatePeaks;
%         save inter
        ga.removeOutlierSubpeaks;
        ga.updateModelAndError;
        
        %calculate integrals of intRanges
        ga.calcIntegrationRanges;
        
%         save inter1
        % swap subpeaks between peaks if needed
        ga.postprocSwapSubpeaks;
        
        ga.correctLocalError(6); 
                
        ga.dispProgress;
        

        %Split overlapping peaks
        ga.postprocPeaksSplittingAlt;
        
        %extract peaks from the tails of peaks
        ga.postprocPeaksTailing;
        
        %adjust overlaping tailing peaks
        ga.postprocPeaksAdjustTails;
        
        %filter small and narrow peaks
        ga.filterPeaks(1.8,12); %1.5
        
        %evaluate peaks for fit quality and overlap
        ga.evalPeaks;
        
        % generate peaks(:).intinfo with integral information about each peak
        ga.updateModelAndError;
        ga.updateIntegrals;
        
        % write general peak info to csv (especially relevant for industry)
        ga.writePeaks;
        
        % write include-range info to csv (especially relevant for industry)
        ga.writeIntRanges;
        
        ga.dispProgress;
        
        if isplot
            %start GUI
            ga.startAnalyzer;
        end
        
        if isREF && ~isknime && isplot
            ga.plotRefComp;
        end
        
        clear ans
        
        
        save([savefol 'gca_res'])
        
        
    case 'analyze'
        % feature for knime
        disp(['gca# Load res-file: ' tempf ])
        load(tempf);
        if exist('ga')
            disp(['gca# ga found in res-file'])
            RESanalysis2(ga);
            ga.dispProgress;
        else
            disp(['gca# ga NOT found in res-file'])
        end
        
     case 'new'
        % feature for knime
        disp(['gca# New analysis'])
        RESanalysis2();
        
    case 'plot'
        % feature for knime
        disp(['gca# Load res-file: ' tempf ])
        load(tempf);
        if exist('ga')
            disp(['gca# ga found in res-file'])
            ga.plot
            ga.dispProgress;
        else
            disp(['gca# ga NOT found in res-file'])
        end
        
        
end

