classdef RESanalysis2 < handle
%% Todo
%


% Properties that correspond to app components
properties (Access = public)
    ver ='0.6.2'

    isdebug=1;

    mainfig         matlab.ui.Figure
    calcfig         matlab.ui.Figure
    mainaxes
    mainaxes_pos
    peaktable

    % tab handles
    hToolGroup

    htemp=[];
    temppeak=[];
    slideval=0;

    tabGroupMain
    tabAdd
    tabMod
    tabRem
    tabDat
    tabFil
    tabInf
    tabs
    chtab

    % ui component handles
    hData
    hAdd
    hMod
    hFil
    hInf

    % axes mode tracker
    axMode='sel'

    listen % legacy variable for detecting view change        
    ifUndoViewRedoView=0; % whether view change caused by undo or redo
    history4UndoView % track last views for undo
    history4RedoView % track last view undo's for redo
    maxHistoryView=20; % = max view undos + 1

    history4Undo % track last actions for undo
    history4Redo % track last undo's for redo
    maxHistory=6; % = max action undos + 1

    multipleClicks = []; % Counts multiple clicks to avoid multiple callbacks
    checkDoubleclick = []; % Detects double clicks for choosing peak in plot

    isSelectCB=1;

    cc =[0    0.4470    0.7410
        0.8500    0.3250    0.0980
        0.9290    0.6940    0.1250
        0.4940    0.1840    0.5560
        0.4660    0.6740    0.1880
        0.3010    0.7450    0.9330
        0.6350    0.0780    0.1840];
end


properties (Access = public)
    ga              GCanalysis
    peakMaxLoc
    hPeak
    hSpec
    hModProf
    currentSel = 0;

end

%% App creation and deletion
methods

    %% CONSTRUCTOR
    function obj = RESanalysis2(ga)
        warning off
        addpath('class')
        addpath('bin')
        addpath('bin_gui')
        addpath(['resource' filesep 'toolstrip_icons'])
        warning on
        disp(['GCA result analyzer - Version ' obj.ver])


        % Create figure and components
        % contains everything that is data !INdependent!
        obj.createComponents;

        % assign GCanalysis object
        if nargin==1
            obj.ga=ga;

            % initialization
            % contains everything that is data dependent
            % after loading a new ga object, just run startup
            obj.startupFcn;
            obj.addToHistory; % Update history for undo
        else
            obj.loadMat;
            if isempty(obj.ga)
               obj.delete 
            end
        end




    end

    %% DESTRUCTOR
    function delete(obj)    
        p = properties(obj);
        for i=1:length(p)
            if isa(eval(['obj.' p{i} ]),'char')
                eval(['obj.' p{i} '=[];'])
            else
                try
                    eval(['delete(obj.' p{i} ')'])
                catch
                    eval(['obj.' p{i} '=[];'])
                end
            end
        end
        delete(obj.mainfig)
        clear obj
    end

    %% main figure close req
    function closeMain_fcn(obj,src,evt)
        delete(obj);
    end
end


methods (Access = private)

    %% Create figure and components
    % Contains everything that is data independent
    function createComponents(obj)
        obj.hToolGroup = matlab.ui.internal.desktop.ToolGroup(['GCanalysis - Peak Postprocessing v' obj.ver]);

        obj.hToolGroup.disableDataBrowser();

        % Store toolgroup reference handle so that app will stay in memory
        jToolGroup = obj.hToolGroup.Peer;
        internal.setJavaCustomData(jToolGroup, obj.hToolGroup);

        if isdeployed
            % jPos: Abstand von linker oberer Ecken
            jPos = com.mathworks.widgets.desk.DTLocation.createExternal(uint16(100), uint16(100), uint16(1200), uint16(170));
            jToolGroup.setPosition(jPos)
        end

        % Create two figures and dock them into the ToolGroup
        obj.mainfig = figure('Name','Peak Analysis','Visible','on');
        obj.mainfig.NumberTitle='off';
        %             obj.hToolGroup.addFigure(obj.mainfig);
        if isdeployed
            screenS=get(0,'ScreenSize');
            sHi=screenS(4);
            set(obj.mainfig,'Position',[100 sHi-20-100-170-700 1200 700])
        end


        obj.hToolGroup.hideViewTab;  % toolstrip View tab is still visible at this point

        import matlab.ui.internal.toolstrip.*  % for convenience below


        obj.tabGroupMain = TabGroup();
        obj.tabDat = Tab('Data');
        obj.tabGroupMain.add(obj.tabDat);
        obj.tabAdd = Tab('Add/Remove');
        obj.tabGroupMain.add(obj.tabAdd);
        obj.tabMod = Tab('Modify');
        obj.tabGroupMain.add(obj.tabMod);
        obj.tabFil = Tab('Filter');
        obj.tabGroupMain.add(obj.tabFil);
        obj.tabInf = Tab('Info');
        obj.tabGroupMain.add(obj.tabInf);

        obj.hToolGroup.addTabGroup(obj.tabGroupMain);


        % add view control and undo to all tabs

        obj.tabs={obj.tabDat,obj.tabMod,obj.tabAdd,obj.tabFil,obj.tabInf};
        obj.chtab = {};

        for i=1:length(obj.tabs)

            hSection = Section('View');
            obj.tabs{i}.add(hSection);

            hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',80);  % add to section as the last column
            obj.chtab{i}.btn_resView = Button('Reset View');
            obj.chtab{i}.btn_resView.Icon=Icon(javax.swing.ImageIcon('Plot_24.png'));
            obj.chtab{i}.btn_resView.ButtonPushedFcn={@(hObj,evnt) obj.resView(hObj,evnt)};
            hColumn.add(obj.chtab{i}.btn_resView);

            hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',80);  % add to section as the last column.

            obj.chtab{i}.btn_undoView = Button('Undo View');
            obj.chtab{i}.btn_undoView.Enabled = false;
            obj.chtab{i}.btn_undoView.Icon=Icon(javax.swing.ImageIcon('Back_16.png'));
            obj.chtab{i}.btn_undoView.ButtonPushedFcn={@(hObj,evnt) obj.undoView(hObj,evnt)};
            hColumn.add(obj.chtab{i}.btn_undoView);

            obj.chtab{i}.btn_redoView = Button('Redo View');
            obj.chtab{i}.btn_redoView.Enabled = false;
            obj.chtab{i}.btn_redoView.Icon=Icon(javax.swing.ImageIcon('Forward_16.png'));
            obj.chtab{i}.btn_redoView.ButtonPushedFcn={@(hObj,evnt) obj.redoView(hObj,evnt)};
            hColumn.add(obj.chtab{i}.btn_redoView);


            hColumn = hSection.addColumn();

            obj.chtab{i}.rad_sel = Button('Select',Icon.SELECT_16);
            obj.chtab{i}.rad_sel.ButtonPushedFcn = {@(hObj,evnt) obj.mouseMode(hObj,evnt,'sel')};
            obj.chtab{i}.rad_sel.Icon=Icon(javax.swing.ImageIcon('Select_16.png'));
            hColumn.add(obj.chtab{i}.rad_sel);

            obj.chtab{i}.rad_pan = Button('Pan',Icon.PAN_16);
            obj.chtab{i}.rad_pan.ButtonPushedFcn = {@(hObj,evnt) obj.mouseMode(hObj,evnt,'pan')};
            obj.chtab{i}.rad_pan.Icon=Icon(javax.swing.ImageIcon('Pan_16.png'));
            hColumn.add(obj.chtab{i}.rad_pan);

            hColumn = hSection.addColumn();

            obj.chtab{i}.rad_zoom = Button('Zoom in',Icon.ZOOM_IN_16);
            obj.chtab{i}.rad_zoom.ButtonPushedFcn = {@(hObj,evnt) obj.mouseMode(hObj,evnt,'zin')};
            obj.chtab{i}.rad_zoom.Icon=Icon(javax.swing.ImageIcon('Zoom_In_16.png'));
            hColumn.add(obj.chtab{i}.rad_zoom);

            obj.chtab{i}.rad_zoomout = Button('Zoom out',Icon.ZOOM_OUT_16);
            obj.chtab{i}.rad_zoomout.ButtonPushedFcn = {@(hObj,evnt) obj.mouseMode(hObj,evnt,'zout')};
            obj.chtab{i}.rad_zoomout.Icon=Icon(javax.swing.ImageIcon('Zoom_Out_16.png'));
            hColumn.add(obj.chtab{i}.rad_zoomout);

            hSection = Section('Undo');
            obj.tabs{i}.add(hSection);

            hColumn = hSection.addColumn();

            obj.chtab{i}.btn_undo = Button('Undo',Icon.UNDO_16);
            obj.chtab{i}.btn_undo.Enabled = false;
            obj.chtab{i}.btn_undo.ButtonPushedFcn = {@(hObj,evnt) obj.undo(hObj,evnt)};
            obj.chtab{i}.btn_undo.Icon=Icon(javax.swing.ImageIcon('Undo_16.png'));
            hColumn.add(obj.chtab{i}.btn_undo);

            obj.chtab{i}.btn_redo = Button('Redo',Icon.REDO_16);
            obj.chtab{i}.btn_redo.Enabled = false;
            obj.chtab{i}.btn_redo.ButtonPushedFcn = {@(hObj,evnt) obj.redo(hObj,evnt)};
            obj.chtab{i}.btn_redo.Icon=Icon(javax.swing.ImageIcon('Redo_16.png'));
            hColumn.add(obj.chtab{i}.btn_redo);


        end


        % INFO tab
        hSection = Section('Debugging');
        obj.tabInf.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',90);  % add to section as the last column
        obj.hData.btn_errorRep = Button('Error Report');
        obj.hData.btn_errorRep.Icon=Icon(javax.swing.ImageIcon('Tools_24.png'));
        obj.hData.btn_errorRep.ButtonPushedFcn={@(hObj,evnt) obj.errorRep(hObj,evnt)};
        hColumn.add(obj.hData.btn_errorRep);


        hSection = Section('Dataset information');
        obj.tabInf.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','right', 'Width',70);  % add to section as the last column
        hColumn.add(Label('File name:'));
        hColumn.add(Label('File folder:'));
        hColumn.add(Label('File path:'));

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',220);  % add to section as the last column
        obj.hInf.lab_name=Label('...');
        hColumn.add(obj.hInf.lab_name);
        obj.hInf.lab_fold=Label('...');
        hColumn.add(obj.hInf.lab_fold);
        obj.hInf.lab_path=EditField('...');
        obj.hInf.lab_path.Editable=false;
        hColumn.add(obj.hInf.lab_path);


        hSection = Section('Loaded result file');
        obj.tabInf.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','right', 'Width',70);  % add to section as the last column
        hColumn.add(Label('File name:'));
        hColumn.add(Label('File folder:'));
        hColumn.add(Label('File path:'));

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',220);  % add to section as the last column
        obj.hInf.lab_nameLoaded=Label('no .mat file loaded');
        hColumn.add(obj.hInf.lab_nameLoaded);
        obj.hInf.lab_foldLoaded=Label('no .mat file loaded');
        hColumn.add(obj.hInf.lab_foldLoaded);
        obj.hInf.lab_pathLoaded=EditField('no .mat file loaded');
        obj.hInf.lab_pathLoaded.Editable=false;
        hColumn.add(obj.hInf.lab_pathLoaded);



        % FILTER tab
        hSection = Section('HEIGHT FILTER');
        obj.tabFil.add(hSection);

        hColumn = hSection.addColumn('Width',100,'HorizontalAlignment','right');
        hColumn.add(Label('Peak Height:'))
        hColumn.add(Label('Affected Peaks:'))

        hColumn = hSection.addColumn('Width',70,'HorizontalAlignment','left');
        obj.hFil.edit_height=EditField('1');
        obj.hFil.edit_height.ValueChangedFcn={@(hobj,evnt) obj.filterPeakHeightEdit(hobj,evnt)};
        hColumn.add(obj.hFil.edit_height);
        obj.hFil.lab_height=Label('...');
        hColumn.add(obj.hFil.lab_height);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',70);  % add to section as the last column
        obj.hFil.btn_remPeak = Button('Remove');
        obj.hFil.btn_remPeak.Icon=Icon(javax.swing.ImageIcon('Delete_24.png'));
        obj.hFil.btn_remPeak.ButtonPushedFcn={@(hobj,evnt) obj.filterPeakHeight(hobj,evnt)};
        hColumn.add(obj.hFil.btn_remPeak);


        hSection = Section('AREA FILTER');
        obj.tabFil.add(hSection);

        hColumn = hSection.addColumn('Width',100,'HorizontalAlignment','right');
        hColumn.add(Label('Peak Area:'))
        hColumn.add(Label('Affected Peaks:'))

        hColumn = hSection.addColumn('Width',70,'HorizontalAlignment','left');
        obj.hFil.edit_int=EditField('1');
        obj.hFil.edit_int.ValueChangedFcn={@(hobj,evnt) obj.filterPeakIntEdit(hobj,evnt)};
        hColumn.add(obj.hFil.edit_int);
        obj.hFil.lab_int=Label('...');
        hColumn.add(obj.hFil.lab_int);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',70);  % add to section as the last column
        obj.hFil.btn_remPeakInt = Button('Remove');
        obj.hFil.btn_remPeakInt.Icon=Icon(javax.swing.ImageIcon('Delete_24.png'));
        obj.hFil.btn_remPeakInt.ButtonPushedFcn={@(hobj,evnt) obj.filterPeakInt(hobj,evnt)};
        hColumn.add(obj.hFil.btn_remPeakInt);


        hSection = Section('Range FILTER');
        obj.tabFil.add(hSection);

%             len=max(obj.ga.modePF.timegrid)-min(obj.ga.modePF.timegrid);
%             timemin=min(obj.ga.modePF.timegrid);
%             timemax=max(obj.ga.modePF.timegrid);

        hColumn = hSection.addColumn('Width',100,'HorizontalAlignment','right');
        hColumn.add(Label('Left Bound:'))
        hColumn.add(Label('Right Bound:'))
        hColumn.add(Label('Affected Peaks:'))

        hColumn = hSection.addColumn('Width',70,'HorizontalAlignment','left');
        % spinner LB
        obj.hFil.spinLB_int = Spinner();
        %obj.hFil.spinLB_int = Spinner([timemin timemax], timemin+0.1*len);  
        obj.hFil.spinLB_int.NumberFormat = 'double';
        obj.hFil.spinLB_int.DecimalFormat = '4f';
%             obj.hFil.spinLB_int.StepSize = (timemax-timemin)/300;
        obj.hFil.spinLB_int.ValueChangedFcn = {@(hobj,evnt) obj.filterPeakRangeEdit(hobj,evnt)};
        hColumn.add(obj.hFil.spinLB_int);
        % spinner RB
        obj.hFil.spinRB_int = Spinner();
%             obj.hFil.spinRB_int = Spinner([timemin timemax], timemax-0.1*len); 
        obj.hFil.spinRB_int.NumberFormat = 'double';
        obj.hFil.spinRB_int.DecimalFormat = '4f';
%             obj.hFil.spinRB_int.StepSize = (timemax-timemin)/300;
%             obj.hFil.spinRB_int.ValueChangedFcn = {@(hobj,evnt) obj.filterPeakRangeEdit(hobj,evnt)};
        hColumn.add(obj.hFil.spinRB_int);
        %label
        obj.hFil.lab_range=Label('...');
        hColumn.add(obj.hFil.lab_range);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',70);  % add to section as the last column
        obj.hFil.btn_remPeakRange = Button('Remove');
        obj.hFil.btn_remPeakRange.Icon=Icon(javax.swing.ImageIcon('Delete_24.png'));
        obj.hFil.btn_remPeakRange.ButtonPushedFcn={@(hobj,evnt) obj.filterPeakRange(hobj,evnt)};
        hColumn.add(obj.hFil.btn_remPeakRange);




        % ADD/REMOVE tab
        hSection = Section('');
        obj.tabAdd.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',100);  % add to section as the last column
        obj.hAdd.btn_addPeak = Button('Add/Draw rect.');
        obj.hAdd.btn_addPeak.Icon=Icon(javax.swing.ImageIcon('Add_24.png'));
        obj.hAdd.btn_addPeak.ButtonPushedFcn={@(hobj,evnt) obj.addPeak(hobj,evnt)};
        hColumn.add(obj.hAdd.btn_addPeak);

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',100);  % add to section as the last column
        obj.hAdd.btn_addPeakFit = Button('Fit');
        obj.hAdd.btn_addPeakFit.Enabled=false;
        obj.hAdd.btn_addPeakFit.Icon=Icon(javax.swing.ImageIcon('Run_16.png'));
        obj.hAdd.btn_addPeakFit.ButtonPushedFcn={@(hobj,evnt) obj.addPeakFit(hobj,evnt)};%%
        hColumn.add(obj.hAdd.btn_addPeakFit);
        obj.hAdd.btn_addPeakOk = Button('Confirm');
        obj.hAdd.btn_addPeakOk.Enabled=false;
        obj.hAdd.btn_addPeakOk.Icon=Icon(javax.swing.ImageIcon('Confirm_16.png'));
        obj.hAdd.btn_addPeakOk.ButtonPushedFcn={@(hobj,evnt) obj.addPeakOk(hobj,evnt)};%%
        hColumn.add(obj.hAdd.btn_addPeakOk);
        obj.hAdd.btn_addPeakCancel = Button('Cancel');
        obj.hAdd.btn_addPeakCancel.Enabled=false;
        obj.hAdd.btn_addPeakCancel.Icon=Icon(javax.swing.ImageIcon('Close_16.png'));
        obj.hAdd.btn_addPeakCancel.ButtonPushedFcn={@(hobj,evnt) obj.addPeakCancel(hobj,evnt)};%%
        hColumn.add(obj.hAdd.btn_addPeakCancel);


        hSection = Section('');
        obj.tabAdd.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',70);  % add to section as the last column
        obj.hAdd.btn_remPeak = Button('Remove');
        obj.hAdd.btn_remPeak.Icon=Icon(javax.swing.ImageIcon('Delete_24.png'));
        obj.hAdd.btn_remPeak.ButtonPushedFcn={@(hobj,evnt) obj.deleteSelPeak(hobj,evnt)};
        hColumn.add(obj.hAdd.btn_remPeak);

        % MODIFY tab
        hSection = Section('Current Peak Info');
        obj.tabMod.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','right', 'Width',100);  % add to section as the last column
        hColumn.add(Label('Peak #:'));
        hColumn.add(Label('Subpeaks:'));

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',50);  % add to section as the last column
        obj.hMod.lab_peakidx= Label('...');
        hColumn.add(obj.hMod.lab_peakidx);

        obj.hMod.lab_submods= Label('...');
        hColumn.add(obj.hMod.lab_submods);


        hSection = Section('Splitting');
        obj.tabMod.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',70);  % add to section as the last column
        obj.hMod.btn_split = Button('Split Peak');
        obj.hMod.btn_split.Icon=Icon(javax.swing.ImageIcon('Cut_24.png'));
        obj.hMod.btn_split.ButtonPushedFcn={@(hobj,evnt) obj.splitPeak(hobj,evnt)};
        hColumn.add(obj.hMod.btn_split);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',200);  % add to section as the last column
        obj.hMod.time_slide=Slider([1 2], 1);
        obj.hMod.time_slide.Enabled=false;
        obj.hMod.time_slide.ValueChangedFcn={@(hobj,evnt) obj.splitChange(hobj,evnt)};
        hColumn.add(obj.hMod.time_slide);

        obj.hMod.btn_splitOk = Button('Confirm');
        obj.hMod.btn_splitOk.Icon=Icon(javax.swing.ImageIcon('Confirm_16.png'));
        obj.hMod.btn_splitOk.ButtonPushedFcn={@(hobj,evnt) obj.splitOk(hobj,evnt)};
        obj.hMod.btn_splitOk.Enabled=false;
        hColumn.add(obj.hMod.btn_splitOk);

        obj.hMod.btn_splitCancel = Button('Cancel');
        obj.hMod.btn_splitCancel.Icon=Icon(javax.swing.ImageIcon('Close_16.png'));
        obj.hMod.btn_splitCancel.ButtonPushedFcn={@(hobj,evnt) obj.splitCancel(hobj,evnt)};
        obj.hMod.btn_splitCancel.Enabled=false;
        hColumn.add(obj.hMod.btn_splitCancel);

        %Refit
        hSection = Section('Refit existing peak');
        obj.tabMod.add(hSection); 

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',100);  
        obj.hMod.btn_PeakReFit = Button('Refit selected');
        obj.hMod.btn_PeakReFit.Enabled=true;
        obj.hMod.btn_PeakReFit.Icon=Icon(javax.swing.ImageIcon('Run_24.png'));
        obj.hMod.btn_PeakReFit.ButtonPushedFcn={@(hobj,evnt) obj.ReFitPeak(hobj,evnt)};%%
        hColumn.add(obj.hMod.btn_PeakReFit);

        hColumn=hSection.addColumn('HorizontalAlignment','left', 'Width',100); 
        obj.hMod.btn_PeakReFitOk = Button('Confirm');
        obj.hMod.btn_PeakReFitOk.Enabled=false;
        obj.hMod.btn_PeakReFitOk.Icon=Icon(javax.swing.ImageIcon('Confirm_16.png'));
        obj.hMod.btn_PeakReFitOk.ButtonPushedFcn={@(hobj,evnt) obj.ReFitOk(hobj,evnt)};%%
        hColumn.add(obj.hMod.btn_PeakReFitOk);
        obj.hMod.btn_PeakReFitCancel = Button('Cancel');
        obj.hMod.btn_PeakReFitCancel.Enabled=false;
        obj.hMod.btn_PeakReFitCancel.Icon=Icon(javax.swing.ImageIcon('Close_16.png'));
        obj.hMod.btn_PeakReFitCancel.ButtonPushedFcn={@(hobj,evnt) obj.ReFitCancel(hobj,evnt)};%%
        hColumn.add(obj.hMod.btn_PeakReFitCancel);



        % DATA tab
        hSection = Section('FILE');
        obj.tabDat.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',50);  % add to section as the last column
        obj.hData.btn_loadM = Button('Load');
        obj.hData.btn_loadM.Icon=Icon(javax.swing.ImageIcon('Open_24.png'));
        obj.hData.btn_loadM.ButtonPushedFcn={@(hObj,evnt) obj.loadMat(hObj,evnt)};
        hColumn.add(obj.hData.btn_loadM);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',50);  % add to section as the last column
        obj.hData.btn_saveM = Button('Save');
        obj.hData.btn_saveM.Icon=Icon(javax.swing.ImageIcon('Save_24.png'));
        obj.hData.btn_saveM.ButtonPushedFcn={@(hObj,evnt) obj.saveMat(hObj,evnt)};
        hColumn.add(obj.hData.btn_saveM);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',50);  % add to section as the last column
        obj.hData.btn_saveE = Button('Export');
        obj.hData.btn_saveE.Icon=Icon(javax.swing.ImageIcon('Compare_24.png'));
        obj.hData.btn_saveE.ButtonPushedFcn={@(hObj,evnt) obj.exportCsv(hObj,evnt)};
        hColumn.add(obj.hData.btn_saveE);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',85);  % add to section as the last column
        obj.hData.btn_saveME = Button('Save/Export');
        obj.hData.btn_saveME.Icon=Icon(javax.swing.ImageIcon('Save_All_24.png'));
        obj.hData.btn_saveME.ButtonPushedFcn={@(hObj,evnt) obj.save_exp(hObj,evnt)};
        hColumn.add(obj.hData.btn_saveME);

        hSection = Section('Evaluation');
        obj.tabDat.add(hSection);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',80);  % add to section as the last column
        obj.hData.btn_showgaplot = Button('Calc. Info');
        obj.hData.btn_showgaplot.Icon=Icon(javax.swing.ImageIcon('Plot_24.png'));
        obj.hData.btn_showgaplot.ButtonPushedFcn={@(hObj,evnt) obj.calcPlot(hObj,evnt)};
        hColumn.add(obj.hData.btn_showgaplot);

        hColumn=hSection.addColumn('HorizontalAlignment','center', 'Width',80);  % add to section as the last column
        obj.hData.btn_showgatime = Button('Calc. Time');
        obj.hData.btn_showgatime.Icon=Icon(javax.swing.ImageIcon('Continue_24.png'));
        obj.hData.btn_showgatime.ButtonPushedFcn={@(hObj,evnt) obj.calcTime(hObj,evnt)};
        hColumn.add(obj.hData.btn_showgatime);






        % Main axes
        obj.mainaxes = axes('parent',obj.mainfig,'units','pixel','Visible','off');
        title(obj.mainaxes, [])
        xlabel(obj.mainaxes, 'Retentionszeit')
        ylabel(obj.mainaxes, 'Intensitaet')
        obj.mainaxes.Position = [430 60 940 540];

        obj.mainaxes_pos=uicontrol('Style','text',...
            'String','...',...
            'Position',[780 490 120 18],...
            'Backgroundcolor','white',...
            'ForegroundColor',0.4*[1 1 1]);

        % Peak Table
        loadJava
        warning off
        obj.peaktable = uix.widget.Table('Parent',obj.mainfig,'Visible','off');
        warning on
        obj.peaktable.ColumnName={'#';'Position'; 'Integral';'Overlap';'Fit'};
        obj.peaktable.Units='pixels';
        obj.peaktable.Position=[30 30 330 570];
        obj.peaktable.CellSelectionCallback={@(hObj,evnt) obj.cb_table(hObj,evnt)};
        obj.peaktable.CellEditCallback={@(hObj,evnt) obj.cb_table(hObj,evnt)};
        obj.peaktable.ColumnFormat={'char','char','char','char','char'};
        obj.peaktable.ColumnEditable=[false false false false false];
        obj.peaktable.ColumnPreferredWidth=[50 50 100 50 50];


        obj.listen.xlim=[];
        obj.listen.ylim=[];

        obj.history4UndoView={};
        obj.history4RedoView={};

        obj.history4Undo={};
        obj.history4Redo={};

        obj.hToolGroup.open();
        obj.hToolGroup.addFigure(obj.mainfig);

        %resize callback
        set(obj.mainfig,'SizeChangedFcn',{@(src,cbd) obj.sizech_fcn(src,cbd)});
        set(obj.mainfig,'WindowButtonMotionFcn',{@(src,cbd) obj.mousemove_fcn(src,cbd)});
%             set(obj.mainfig,'WindowKeyReleaseFcn',{@(src,cbd) obj.axBtnUp_fcn(src,cbd)});
        addlistener(obj.mainaxes.XRuler,'MarkedClean',@(~,~) obj.axBtnUp_fcn);
%             addlistener(obj.mainaxes,'XLim','PostSet',@(~,~) obj.axBtnUp_fcn);
        set(obj.mainfig,'CloseRequestFcn',{@(src,evt) obj.closeMain_fcn(src,evt)});
        figure(obj.mainfig)

        
    end

    %% Detect mouse release in axis, add view hist entry
    function axBtnUp_fcn(obj,init)
        % Mostly legacy code to detect, when a change of xlim for view
        % history
        debug=0;

        if nargin==1
            init=0;
        end

        if debug
            disp(obj.listen.xlim)
            disp(get(obj.mainaxes,'XLim'))
            disp(' ')
        end

        cxlim=get(obj.mainaxes,'XLim');
        cylim=get(obj.mainaxes,'YLim');
        % init just puts the current limits into limits.xlim/ylim !only!
        if init
            obj.listen.xlim=cxlim;
            obj.listen.ylim=cylim;
            % obj.addToHistoryView(); % ! Don't add to history
            return
        end

        % check for empty listen.xlim while initialization
        if isempty(obj.listen.xlim)
            return
        end

        % check if pan mode or no xlim change
        if isequal(obj.listen.xlim(end,:),cxlim) || strcmp(obj.axMode,'pan') 
            return
        end

        if obj.ifUndoViewRedoView
            obj.ifUndoViewRedoView = 0;
            return
        end

        % legacy: add to limit history
        obj.listen.xlim(end+1,:)=cxlim;
        obj.listen.ylim(end+1,:)=cylim;
        if size(obj.listen.xlim,1)>3
            obj.listen.xlim(1,:)=[];
            obj.listen.ylim(1,:)=[];
        end

        % add to view history
        obj.addToHistoryView();

        if debug
            disp('Saved View')
        end

    end

    %% Callback for end of pan to add the view to view history
    function axBtnUp_endOfPan(obj,init)
        obj.addToHistoryView();
    end


    %% Add current view to view history for undo view
    function addToHistoryView(obj)
        % Save the view position

        % Add current state to undo view stack
        obj.history4UndoView{end+1} = [get(obj.mainaxes,'XLim');get(obj.mainaxes,'YLim')];
        if length(obj.history4UndoView)>obj.maxHistoryView
            obj.history4UndoView = obj.history4UndoView(2:end);
        end
        % Clear redo view stack
        obj.history4RedoView = {};
        % Disable redo view button
        for i=1:length(obj.tabs)
            obj.chtab{i}.btn_redoView.Enabled = false;
        end
        % Enable undo view button if undo view stack contains at least 2 elements
        if length(obj.history4UndoView)>=2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undoView.Enabled = true;
            end
        end
    end

    %% Undo view
    function undoView(obj,hObj,evt)        
        % Last state is the current state
        % Check whether undo view stack contains less than 2 elements
        if length(obj.history4UndoView)<2
            error('Undo View stack empty');
            return
        end
        % Add current state to redo view stack
        obj.history4RedoView{end+1} = obj.history4UndoView{end}; 
        % Enable redo view button
        for i=1:length(obj.tabs)
            obj.chtab{i}.btn_redoView.Enabled = true;
        end
        % Actual undo view
        set(obj.mainaxes,'xlim',obj.history4UndoView{end-1}(1,:),'ylim',obj.history4UndoView{end-1}(2,:))
        % Don't update view history after undo view
        obj.ifUndoViewRedoView = 1;
        % Remove the last state from undo view stack
        obj.history4UndoView = obj.history4UndoView(1:end-1);
        % Disable undo view button if undo view stack contains less than 2 elements
        if length(obj.history4UndoView)<2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undoView.Enabled = false;
            end
        end
    end

    %% Redo view
    function redoView(obj,hObj,evt)
        % Check whether redo view stack is empty
        if isempty(obj.history4RedoView)
            error('Redo View stack empty');
            return
        end
        % Add current state to undo view stack and check length
        obj.history4UndoView{end+1} = obj.history4RedoView{end};
        if length(obj.history4UndoView)>obj.maxHistoryView
            obj.history4UndoView = obj.history4UndoView(2:end);
        end
        % Actual redo view 
        set(obj.mainaxes,'xlim',obj.history4RedoView{end}(1,:),'ylim',obj.history4RedoView{end}(2,:))
        % Don't update view history after redo view
        obj.ifUndoViewRedoView = 1;
        % Remove the last state from redo view stack
        obj.history4RedoView = obj.history4RedoView(1:end-1);
        % Disable redo view button if redo view stack is empty
        if isempty(obj.history4RedoView)
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_redoView.Enabled = false;
            end
        end
        % Enable undo view button if undo view stack contains at least 2 elements
        if length(obj.history4UndoView)>=2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undoView.Enabled = true;
            end
        end
    end



    %% Add current state to history for undo
    function addToHistory(obj)
        % Save the state of the object at the beginning and after each
        % startupFcn (except in the undo and redo functions themself)
        % Add the state of the object to the history after(!) each action
        % The state of the object is saved using an undocumented function
        % to create a deep copy of each handle object.

        % Add current state to undo stack
        obj.history4Undo{end+1} = getByteStreamFromArray(obj.ga);
        if length(obj.history4Undo)>obj.maxHistory
            obj.history4Undo = obj.history4Undo(2:end);
        end
        % Clear redo stack
        obj.history4Redo = {};
        % Disable redo button
        for i=1:length(obj.tabs)
            obj.chtab{i}.btn_redo.Enabled = false;
        end
        % Enable undo button if undo stack contains at least 2 elements
        if length(obj.history4Undo)>=2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undo.Enabled = true;
            end
        end
    end

    %% Undo
    function undo(obj,hObj,evt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        % Last state is the current state
        % Check whether undo stack contains less than 2 elements
        if length(obj.history4Undo)<2
            error('Undo stack empty');
            return
        end
        % Add current state to redo stack
        obj.history4Redo{end+1} = obj.history4Undo{end}; %= getByteStreamFromArray(obj.ga);
        % Enable redo button
        for i=1:length(obj.tabs)
            obj.chtab{i}.btn_redo.Enabled = true;
        end
        % Actual undo
        delete(obj.ga)
        obj.ga = getArrayFromByteStream(obj.history4Undo{end-1});
        obj.startupFcn;
        obj.deselectPeak;
        % Remove the last state from undo stack
        obj.history4Undo = obj.history4Undo(1:end-1);
        % Disable undo button if undo stack contains less than 2 elements
        if length(obj.history4Undo)<2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undo.Enabled = false;
            end
        end
        % Function ends here
        obj.enableButtons();
    end


    %% Redo
    function redo(obj,hObj,evt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        % Check whether redo stack is empty
        if isempty(obj.history4Redo)
            error('Redo stack empty');
            return
        end
        % Add current state to undo stack and check length
        obj.history4Undo{end+1} = obj.history4Redo{end}; % = getByteStreamFromArray(obj.ga);
        if length(obj.history4Undo)>obj.maxHistory
            obj.history4Undo = obj.history4Undo(2:end);
        end
        % Actual redo
        delete(obj.ga);
        obj.ga=getArrayFromByteStream(obj.history4Redo{end});
        obj.startupFcn;
        obj.deselectPeak;
        % Remove the last state from redo stack
        obj.history4Redo = obj.history4Redo(1:end-1);
        % Disable redo button if redo stack is empty
        if isempty(obj.history4Redo)
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_redo.Enabled = false;
            end
        end
        % Enable undo button if undo stack contains at least 2 elements
        if length(obj.history4Undo)>=2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undo.Enabled = true;
            end
        end
        % Function ends here
        obj.enableButtons();
    end


    %% Code that executes after component creation
    % contains everything that is data dependent
    function startupFcn(obj)
        if isempty(obj.ga)
            error('no ga object loaded');
            return
        end
        hdia= msgbox({'Wait for GUI to complete initialization.','This box closes automatically.'},'Info');

        obj.mainfig.Visible='off';
        obj.mainaxes.Visible='off';
        obj.peaktable.set('Visible','off');


        %init spinner ranges
        len=max(obj.ga.modePF.timegrid)-min(obj.ga.modePF.timegrid);
        timemin=min(obj.ga.modePF.timegrid);
        timemax=max(obj.ga.modePF.timegrid);

        %spinnerLB
        obj.hFil.spinLB_int.Limits=[timemin timemax];
        obj.hFil.spinLB_int.Value=timemin+0.1*len;
        obj.hFil.spinLB_int.StepSize = (timemax-timemin)/300;
        obj.hFil.spinLB_int.ValueChangedFcn = {@(hobj,evnt) obj.filterPeakRangeEdit(hobj,evnt)};

        %spinnerRB
        obj.hFil.spinRB_int.Limits=[timemin timemax];
        obj.hFil.spinRB_int.Value=timemax-0.1*len;
        obj.hFil.spinRB_int.StepSize = (timemax-timemin)/300;
        obj.hFil.spinRB_int.ValueChangedFcn = {@(hobj,evnt) obj.filterPeakRangeEdit(hobj,evnt)};

        % Add file name and file path
        if isfield(obj.ga.para,'fName') % Check whether field 'name' exists for backwards compability
            fNameSplit = regexp(obj.ga.para.fName,filesep,'split');
            obj.hInf.lab_name.Text=fNameSplit{end};
            obj.hInf.lab_fold.Text=fNameSplit{end-1};
            obj.hInf.lab_path.Value=obj.ga.para.fName;
        else
            obj.hInf.lab_name.Text='unknown';
            obj.hInf.lab_fold.Text='unknown';
            obj.hInf.lab_path.Value='unknown';
        end

        %ensure backwards compatibility if peaks.intinfo not yet present 
        obj.ga.updateIntegrals;

        %update labels
        obj.filterPeakHeightEdit([],[]);
        obj.filterPeakIntEdit([],[]);
        obj.filterPeakRangeEdit([],[]);

        obj.hSpec=[];
        delete(obj.hModProf);
        obj.hModProf=[];

        delete(obj.hPeak);
        obj.hPeak=[];

        obj.updateTable;
        obj.updatePlot;

        obj.sizech_fcn;
        figure(obj.mainfig)
        obj.mainfig.Visible='on';
        obj.mainaxes.Visible='on';
        obj.peaktable.set('Visible','on');

        obj.peakMaxLoc=zeros(length(obj.ga.peaks),1);
        for i=1:length(obj.ga.peaks)
            obj.peakMaxLoc(i)=obj.ga.peaks(i).peakstruct.maxlocation;
        end

        %click axes callback
        actPlotClick(obj)

        obj.mouseMode([],[],'sel'); % Switch to 'sel', to avoid problems with saving view history

        %init listen.xlim/ylim
        obj.axBtnUp_fcn(1)

        if ishandle(hdia)
            delete(hdia);
        end

    end

    %% CLICKABLE AXIS - NO
    function deactPlotClick(obj)
        set(obj.mainaxes,'ButtonDownFcn',[]);
        set(obj.hPeak,'ButtonDownFcn',[]);
    end

    %% CLICKABLE AXIS - YES
    function actPlotClick(obj)
        set(obj.mainaxes,'ButtonDownFcn',{@(src,cbd) obj.clickAxes(src,cbd)});
        set(obj.hPeak,'ButtonDownFcn',{@(src,cbd) obj.clickAxes(src,cbd)});
    end

    %% CLEAR HTEMP
    function clearHTemp(obj)
        for i=1:length(obj.htemp)
            delete(obj.htemp(i));
        end
        obj.htemp=[];
    end
    
    %% Disable buttons while calculating
    function locked = disableButtons(obj)
        if isempty(obj.multipleClicks) % Single click
            obj.multipleClicks = 1;
            locked = 0;
            % Disable all relevant buttons
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undo.Enabled = false;
                obj.chtab{i}.btn_redo.Enabled = false;
            end    
            obj.hFil.btn_remPeakRange.Enabled=false;
            obj.hFil.btn_remPeakInt.Enabled=false;
            obj.hFil.btn_remPeak.Enabled=false;
            obj.hMod.btn_split.Enabled=false;
            obj.hMod.btn_PeakReFit.Enabled=false;
            obj.hAdd.btn_addPeak.Enabled=false;
            obj.hAdd.btn_remPeak.Enabled=false;
            obj.hData.btn_loadM.Enabled=false;
            obj.hData.btn_saveM.Enabled=false;
            obj.hData.btn_saveE.Enabled=false;
            obj.hData.btn_saveME.Enabled=false;
        else % Multiple click during 1 second
            %multipleClicks = [];
            %pause(1)
            locked = 1;
        end
    end
    
    %% Enable buttons after calculation
    function enableButtons(obj)
        obj.multipleClicks = [];
        % Enable all relevant buttons
        % Enable undo button if undo stack contains at least 2 elements
        if length(obj.history4Undo)>=2
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_undo.Enabled = true;
            end
        end
        % Enable redo button if redo stack is not empty
        if ~isempty(obj.history4Redo)
            for i=1:length(obj.tabs)
                obj.chtab{i}.btn_redo.Enabled = true;
            end
        end
        obj.hFil.btn_remPeakRange.Enabled=true;
        obj.hFil.btn_remPeakInt.Enabled=true;
        obj.hFil.btn_remPeak.Enabled=true;
        obj.hMod.btn_split.Enabled=true;
        obj.hMod.btn_PeakReFit.Enabled=true;
        obj.hAdd.btn_addPeak.Enabled=true;
        obj.hAdd.btn_remPeak.Enabled=true;
        obj.hData.btn_loadM.Enabled=true;
        obj.hData.btn_saveM.Enabled=true;
        obj.hData.btn_saveE.Enabled=true;
        obj.hData.btn_saveME.Enabled=true;
    end
    
    %% Deselect peaks (e.g. after deleting)
    function deselectPeak(obj)
        obj.peaktable.SelectedRows = [];
        obj.currentSel = 0;
    end

    %% SPLIT PEAK
    function splitPeak(obj,hObj,evnt)
        if obj.currentSel==0
            errordlg('Select a peak first.')
            return
        end

        p=obj.ga.peaks(obj.currentSel);
        if size(p.peakstruct.paramat,2)>1
            % Disable buttons
            isLocked = obj.disableButtons();
            if isLocked
                return
            end
            
            set(obj.peaktable,'Enable','off');
            obj.hMod.time_slide.Enabled=true;
            obj.hMod.btn_splitOk.Enabled=true;
            obj.hMod.btn_splitCancel.Enabled=true;

            obj.hMod.time_slide.Limits=[1 size(p.peakstruct.paramat,2)-1];
            obj.hMod.time_slide.Ticks=size(p.peakstruct.paramat,2)-1;

            obj.deactPlotClick
            obj.splitChange([],[],1);

        else
            errordlg('Currently selected peak has only one modelled subpeak. Abort splitting.')
        end

    end

    %% SPLIT PEAK SLIDER CHANGE
    function splitChange(obj,hObj,evnt,init)
        if nargin<=3
            init=0;
        end

        disp('split peak - slideChange: start')
        val=obj.hMod.time_slide.Value;
        if val==obj.slideval && ~init 
            disp('split peak - slideChange: same value interrupt')
            return
        end

        obj.slideval=val;
        p=obj.ga.peaks(obj.currentSel);

        valx=mean(p.peakstruct.paramat(1,[val val+1]));
        ma_global=max(obj.ga.profPF.signal);

        [pl,pr]=p.split(valx,obj.ga.rawP.timegrid,ma_global,0);

        obj.clearHTemp

        obj.htemp(1)=pl.fill('r',obj.mainaxes);
        obj.htemp(2)=pr.fill('b',obj.mainaxes);


    end

    %% SPLIT PEAK CONFIRM
    function splitOk(obj,hObj,evnt)
        set(obj.peaktable,'Enable','on');
        obj.hMod.time_slide.Enabled=false;
        obj.hMod.btn_splitOk.Enabled=false;
        obj.hMod.btn_splitCancel.Enabled=false;
        obj.enableButtons();
        obj.clearHTemp

        %split prepare
        p=obj.ga.peaks(obj.currentSel);
        val=obj.hMod.time_slide.Value;
        valx=mean(p.peakstruct.paramat(1,[val val+1]));
        ma_global=max(obj.ga.profPF.signal);

        %split
        [pl,pr]=p.split(valx,obj.ga.rawP.timegrid,ma_global,0);

        %delete old peak
        cr=obj.currentSel;
        obj.ga.peaks(cr)=[];
        delete(obj.hPeak(cr))   %remove plot
        obj.hPeak(cr)=[];       %delete handle

        % add left and right part
        obj.ga.peaks=[obj.ga.peaks(1:cr-1);pl;pr; obj.ga.peaks(cr:end)];
        obj.ga.evalPeaks(max(1,cr-1):min(cr+2,length(obj.ga.peaks)));

        %update integrals
        obj.ga.updateIntegrals;

        %update pl and pr (fit, overlap, intinfo)
        pl=obj.ga.peaks(cr);
        pr=obj.ga.peaks(cr+1);

        hpl=plot(pl,'w',2,obj.mainaxes);
        hpr=plot(pr,'w',2,obj.mainaxes);
        obj.hPeak=[obj.hPeak(1:cr-1); hpl; hpr; obj.hPeak(cr:end)];
        obj.updateColor;

        %activate Click CBs in Plot
        obj.actPlotClick

        % update table
        obj.peaktable.set('Visible','off');
        obj.peaktable.removeRow(cr);
        obj.peaktable.addRow({num2str(cr,'%i'),num2str(pl.peakstruct.maxlocation,'%.3f'),num2str(pl.intinfo.abs,'%.3e'), num2str(pl.overlap,'%.3f'), num2str(pl.fit,'%.3f')},cr);
        obj.peaktable.addRow({num2str(cr+1,'%i'),num2str(pr.peakstruct.maxlocation,'%.3f'),num2str(pr.intinfo.abs,'%.3e'), num2str(pr.overlap,'%.3f'), num2str(pr.fit,'%.3f')},cr+1);
        for i=cr+2:length(obj.ga.peaks)
            obj.peaktable.setCell(i,1,num2str(i,'%i'));
        end
        obj.peaktable.set('Visible','on');

        %update peak locations
        obj.peakMaxLoc(cr)=[];
        obj.peakMaxLoc=[obj.peakMaxLoc(1:cr-1); pl.peakstruct.maxlocation; pr.peakstruct.maxlocation; obj.peakMaxLoc(cr:end)];

        obj.peaktable.setSelectedRow(cr); % Select correct peak in the table
        obj.addToHistory; % Update history for undo
    end

    %% SPLIT PEAK CANCEL
    function splitCancel(obj,hObj,evnt)
        set(obj.peaktable,'Enable','on');
        obj.hMod.time_slide.Enabled=false;
        obj.hMod.btn_splitOk.Enabled=false;
        obj.hMod.btn_splitCancel.Enabled=false;
        obj.enableButtons();
        obj.clearHTemp
        obj.actPlotClick
    end

    %% SHOW BASELINE/CLUSTER/FIT INFO
    function calcPlot(obj,hObj,evnt)
        obj.calcfig = figure('Name','Calculation results','Visible','off');
        obj.calcfig.set('Position',[100 100 1200 900]);
        obj.calcfig.NumberTitle='off';
        obj.ga.plot(obj.calcfig)
        obj.calcfig.set('Visible','on');
        figure(obj.calcfig)
    end

    %% SHOW CALCULATION TIME
    function calcTime(obj,hObj,evnt)
        c=obj.ga.dispProgress;
        t = cell2table(c,'VariableNames',{'Task','Time'});

        fig = uifigure('Position',[300,600,240,190]);

        uit = uitable(fig,'Position',[10,10,220,170]);
        uit.Data = t;
        uit.ColumnSortable = false;
        uit.ColumnWidth = {'auto',80};
        uit.ColumnEditable = false;

    end

    %% CLICK AXES FCN
    function clickAxes(obj,hObj,evnt)
        % Only catch double clicks on peaks
        % persistent checkDoubleclick % This is alternative to obj.
        if isempty(obj.checkDoubleclick)
            obj.checkDoubleclick = 1;
            pause(0.5); % Delay for double click
            if obj.checkDoubleclick == 1 % Single click
                obj.checkDoubleclick = [];
                % Do nothing
            end
        else % Double click
            obj.checkDoubleclick = [];
            % Choose peak
            loc=get(obj.mainaxes, 'Currentpoint');
            xs=loc(1,1);
            [~,id]=min(abs(obj.peakMaxLoc-xs));
            obj.peaktable.setSelectedRow(id);
            %obj.peaktable.SelectedRows = id;
            
        end
    end

    %% RESET AXES LIMITS
    function resView(obj,hObj,evnt)
        axis(obj.mainaxes,'tight')
        yy=ylim(obj.mainaxes);
        ylim(obj.mainaxes,[min(yy(2),0) yy(2)])
    end

    %% MOUSE MODE SELECT/ZOOM/PAN
    function mouseMode(obj,hObj,evnt,type)
        switch type
            case 'sel'
                obj.axMode='sel';
                zoom off
                pan off
                set(obj.mainaxes,'ButtonDownFcn',{@(src,cbd) obj.clickAxes(src,cbd)});
            case 'pan'
                obj.axMode='pan';
                %pan(obj.mainaxes)
                h=pan(obj.mainfig);
                h.Enable='on';
                h.ActionPostCallback = @(~,~) obj.axBtnUp_endOfPan;
            case 'zin'
                obj.axMode='zin';
                h=zoom(obj.mainaxes);
                h.Enable='on';
                h.Direction='in';
            case 'zout'
                obj.axMode='zout';
                h=zoom(obj.mainaxes);
                h.Enable='on';
                h.Direction='out';

        end
        %
        %             if obj.hData.rad_sel.Value
        %
        %             elseif obj.hData.rad_zoom.Value
        %                 h=zoom(obj.mainaxes);
        %                 h.Enable='on';
        %                 h.Direction='in';
        %             elseif obj.hData.rad_zoomout.Value
        %                 h=zoom(obj.mainaxes);
        %                 h.Enable='on';
        %                 h.Direction='out';
        %             elseif obj.hData.rad_pan.Value
        %                 pan(obj.mainaxes)
        %             end
    end

    %% update graphical elements
    function updatePlot(obj)            
        obj.hSpec=plot(obj.ga.profPF,'k',1,obj.mainaxes);
        hold(obj.mainaxes,'on');
        obj.hModProf=plot(obj.ga.modePF,'b',1,obj.mainaxes);
        axis(obj.mainaxes,'tight')
        yy=ylim(obj.mainaxes);
        ylim(obj.mainaxes,[min(yy(2),0) yy(2)])

        obj.hPeak=zeros(length(obj.ga.peaks),1);
        id=1;
        for i=1:length(obj.ga.peaks)
            obj.hPeak(i)=plot(obj.ga.peaks(i),obj.cc(id,:),2,obj.mainaxes);
            id=id+1;
            if id==8
                id=1;
            end
        end
        grid(obj.mainaxes, 'on')
    end

    function updateColor(obj)
        id=1;
        for i=1:length(obj.hPeak)
            set(obj.hPeak(i),'Color',obj.cc(id,:));
            id=id+1;
            if id==8
                id=1;
            end
        end
    end

    %% UPDATE TABLE
    function updateTable(obj)
        dataCell=cell(length(obj.ga.peaks),5);
        for i=1:length(obj.ga.peaks)
            dataCell{i,1}=num2str(i,'%i');
            dataCell{i,2}=num2str(obj.ga.peaks(i).peakstruct.maxlocation,'%.3f');
            dataCell{i,3}=num2str(obj.ga.peaks(i).intinfo.abs,'%.3e');
            dataCell{i,4}=num2str(obj.ga.peaks(i).overlap,'%.3f');
            dataCell{i,5}=num2str(obj.ga.peaks(i).fit,'%.3f');
        end
        obj.peaktable.set('Visible','off');
        obj.peaktable.Data=dataCell;
        obj.peaktable.set('Visible','on');
    end
end

methods (Access = private)
    %% LOAD MAT
    function obj = loadMat(obj,hObj,evnt)
        if ~isempty(obj.ga)
            answer = questdlg('Proceed? Unsaved results will be lost!', ...
                'Warning', ...
                'Proceed','Abort','Abort');
            % Handle response
            if strcmpi(answer,'Abort')
                return
            end
        end

        [file,path] = uigetfile('*.mat','Select a GCD result file');
        if isequal(file,0)
            return
        end

        cload=load(fullfile(path,file),'-mat');
        if isfield(cload,'ga')
            %delete old ga
            delete(obj.ga);
            obj.ga=cload.ga;

            % Add file name and file path
            fNameSplit = regexp(fullfile(path,file),filesep,'split');
            obj.hInf.lab_nameLoaded.Text=fNameSplit{end};
            obj.hInf.lab_foldLoaded.Text=fNameSplit{end-1};
            obj.hInf.lab_pathLoaded.Value=fullfile(path,file);
        else
            warndlg('Selected file contains no variable "ga".')
            return
        end


        obj.startupFcn;
        obj.addToHistory; % Update history for undo
    end

    %% SAVE MAT
    function obj = saveMat(obj,hObj,evnt)
        [fname,fpath] = uiputfile('*.mat','Select file');
        if max(fname)>0
            ga=obj.ga;
            save([fpath fname],'ga');
        end
    end

    %% EXPORT CSV
    function obj = exportCsv(obj,hObj,evnt)
        [fname,fpath] = uiputfile('*.csv','Select file');
        if max(fname)>0
            obj.ga.writePeaks([fpath fname]);
        end
    end

    %% SAVE & EXPORT
    function obj = save_exp(obj,hObj,evnt)
        [fname,fpath] = uiputfile('*.mat;*.csv','Select file');
        if max(fname)>0
            if strcmpi(fname(end-3:end),'.csv')
                fname=fname(end-2:end);
            end
            if strcmpi(fname(end-3:end),'.mat')
                fname=fname(end-2:end);
            end
            ga=obj.ga;
            save([fpath fname '.mat'],'ga');

            obj.ga.writePeaks([fpath fname '.csv']);
        end
    end

    %% FILTER VIA INT EDIT CHANGE
    function obj = filterPeakIntEdit(obj,hObj,evnt)
        delidx=getDelidx(obj,'integral');
        obj.hFil.lab_int.Text=num2str(sum(delidx));
    end

    %% FILTER VIA INT
    function obj = filterPeakInt(obj,hObj,evnt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        delidx=getDelidx(obj,'integral');
        obj.deletePeakRange(find(delidx));
        obj.deselectPeak;
        % Update history for undo
        obj.addToHistory;
        %update affected peak label
        obj.filterPeakIntEdit([],[]);
        % Function ends here
        obj.enableButtons();
    end

    %% FILTER VIA HEIGHT EDIT CHANGE
    function obj = filterPeakHeightEdit(obj,hObj,evnt)
        delidx=getDelidx(obj,'height');
        obj.hFil.lab_height.Text=num2str(sum(delidx));
    end

    %% AUX - get delete indices for filter
    function delidx=getDelidx(obj,type)
        switch type
            case 'height'
                % get edit height val
                sval=obj.hFil.edit_height.Value;
            case 'integral'
                % get edit height val
                sval=obj.hFil.edit_int.Value;
        end
        val=str2double(sval);
        if isnan(val)
            error('Value in input field is not numeric')
        end

        if val<0
            error('Value in input field is <0, but has to be >=0')
        end

        p=obj.ga.peaks;
        delidx=zeros(size(p));
        for i=1:length(p)
            switch type
                case 'height'
                    if p(i).peakstruct.height <= val
                        delidx(i)=1;
                    end
                case 'integral'
                    if p(i).intinfo.abs <= val
                        delidx(i)=1;
                    end
            end
        end
    end

    %% FILTER VIA HEIGHT
    function obj = filterPeakHeight(obj,hObj,evnt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        delidx=getDelidx(obj,'height');
        obj.deletePeakRange(find(delidx));
        obj.deselectPeak;
        % Update history for undo
        obj.addToHistory;
        %update affected peak label
        obj.filterPeakHeightEdit([],[]);
        % Function ends here
        obj.enableButtons();
    end

    %% AUX - get delete indices for filter
    function delidx=getDelidxRange(obj,lb,rb)
        p=obj.ga.peaks;
        delidx=zeros(size(p));
        for i=1:length(p)
            if p(i).peakstruct.maxlocation <= rb && p(i).peakstruct.maxlocation >= lb
                delidx(i)=1;
            end
        end
    end

    %% FILTER VIA RANGE EDIT CHANGE
    function obj = filterPeakRangeEdit(obj,hObj,evnt)
        lb=obj.hFil.spinLB_int.Value;
        rb=obj.hFil.spinRB_int.Value;
        if lb>rb
            lb=rb;
            obj.hFil.spinLB_int.Value=rb;
        else
            timemin=min(obj.ga.modePF.timegrid);
            timemax=max(obj.ga.modePF.timegrid);
            obj.hFil.spinLB_int.Limits=[timemin rb];
            obj.hFil.spinRB_int.Limits=[lb timemax];
        end

        delidx=getDelidxRange(obj,lb,rb);
        obj.hFil.lab_range.Text=num2str(sum(delidx));
    end

    %% FILTER VIA RANGE
    function obj = filterPeakRange(obj,hObj,evnt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        lb=obj.hFil.spinLB_int.Value;
        rb=obj.hFil.spinRB_int.Value;
        delidx=getDelidxRange(obj,lb,rb);
        obj.deletePeakRange(find(delidx));
        obj.deselectPeak;
        % Update history for undo
        obj.addToHistory;
        %update affected peak label
        obj.filterPeakRangeEdit([],[]);
        % Function ends here
        obj.enableButtons();
    end

    %% DELETE PEAK (Range)
    function obj = deletePeakRange(obj,range)

        %eliminate multiples
        range=unique(range);

        % check input
        validateattributes(range, {'double','int32'},{'integer','>=',1,'<=',length(obj.ga.peaks)})

        mod=zeros(size(obj.ga.modePF.timegrid));
        for i=1:length(range)
            mod=mod+obj.ga.peaks(range(i)).getmod(obj.ga.modePF.timegrid);
        end
        obj.ga.modePF.signal=obj.ga.modePF.signal-mod;

        %delete peak plot
        obj.ga.peaks(range)=[];
        obj.peakMaxLoc(range)=[];
        delete(obj.hPeak(range));
        obj.hPeak(range)=[];

        %update integrals
        obj.ga.updateIntegrals;

        %update model
        delete(obj.hModProf);
        obj.hModProf=plot(obj.ga.modePF,'b',1,obj.mainaxes);

        %update table
        obj.peaktable.set('Visible','off');
        obj.isSelectCB=0;
        obj.updateTable
        obj.currentSel=min(obj.currentSel,length(obj.hPeak));
        ev.Indices=[obj.currentSel 1];
        obj.cb_table([],ev)
        obj.peaktable.set('Visible','on');
    end

    %% DELETE PEAK
    function obj = deleteSelPeak(obj,hObj,evnt)
        if obj.currentSel==0
            errordlg('Select a peak first.')
            return
        end
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        obj.deletePeakRange(obj.currentSel);
        obj.deselectPeak;
        % Update history for undo
        obj.addToHistory;
        % Function ends here
        obj.enableButtons();
    end
    
    %% ADD PEAK
    function obj = addPeak(obj,hObj,evnt)
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        obj.htemp(end+1) = drawrectangle();
        obj.hAdd.btn_addPeakFit.Enabled=true;
        obj.hAdd.btn_addPeakOk.Enabled=true;
        obj.hAdd.btn_addPeakCancel.Enabled=true;
    end

    %% ADD PEAK Fit
    function obj = addPeakFit(obj,hObj,evnt)
        if ~isempty(obj.htemp)
            if ishandle(obj.htemp(1))
                Position=get(obj.htemp(1),'Position');
                lb=Position(1);
                rb=Position(1)+Position(3);


                gx=obj.ga.profPF.timegrid;
                gc=obj.ga.profPF.signal;

                lb=max(lb,gx(1));
                rb=min(rb,gx(end));
                if lb>rb
                    warning('add peak - fit: swap border of roi')
                end

                %find indices
                lbx=find(gx<=lb,1,'last');
                rbx=find(gx>=rb,1,'first');

                if abs(lbx-rbx)<10
                    warning('add peak - fit: empty/too small roi (min of 10 data points)')
                    msgbox('Warning - rectangle is too small (min of 10 data points)')
                    return
                end

                %get gc and model window
                gxt=gx(lbx:rbx);
                gct=gc(lbx:rbx);
                modt=obj.ga.modePF.signal(lbx:rbx);

                rest=gct-modt;

                %Scaling to max 1 / get center idx
                [scal1,idce]=max(rest);
                rest_scal=rest/scal1;
                center=gxt(idce);

                %center to max position
                gxt_shift=gxt-center;

                % approximated width
                widthL=0.15*(rb-lb);
                widthR=0.15*(rb-lb);

                hiscale=0.9;

                %full active window
                actx=false(size(gxt_shift));
                help=linspace(0,1,length(actx));
                win_width=0.2;
                actx(help>=win_width)=true;
                actx(help>(1-win_width))=false;
                actx_w=true(size(gxt_shift));

                if obj.isdebug
                    ft=figure;
                    subplot(2,2,1)
                    plot(gxt,gct,'k',gxt, modt,'b','linewidth',2);
                    legend('signal','model');
                    yy=ylim;

                    subplot(2,2,2)
                    plot(gxt_shift,rest_scal,'k','linewidth',2);
                    hold on
                    plot(gxt_shift([1 end]),obj.ga.rescl.noiselvl/scal1*[1 1],'r--','linewidth',2)
                    plot(0,hiscale,'g*','linewidth',2);
                    plot([-widthL widthR],0.5*[1 1],'g-*','linewidth',2);
                    plot(gxt_shift(actx),zeros(size(gxt_shift(actx))),'m-','linewidth',2)
                    legend('residual','noise level','max','width','data fit win');
                    title('shifted/scaled')
                end


                % set initial values and lower/upper bounds
                p0=[0,...           center
                    1.3*widthL,...   sigma left  (1.3*)
                    1.3*widthR,...   sigma right (1.3*)
                    hiscale...         height
                    0.3 ...         linearity left
                    0.3 ...         linearity right
                    ];
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
                we=[1 10];
                fun2=@(p) tar_opt(p,gxt_shift,rest_scal,1,6,we,[],[],actx,actx_w,0);
                options=optimoptions(@lsqnonlin,...
                    'Display','off',...
                    'MaxIter',100,...
                    'MaxFunEvals',10000,...
                    'TolFun',1e-10,...
                    'TolX',1e-10);
                pOpt=lsqnonlin(fun2,p0,pL,pR,options);

                %model eval / rescale / reshift
                mod_gxt_shift=mod_ga(gxt_shift,pOpt(:));
                pOpt(4)=scal1*pOpt(4);
                pOpt(1)=pOpt(1)+center;

                fullpeak=mod_ga(gx,pOpt(:));

                if obj.isdebug
                    figure(ft);
                    subplot(2,2,3)
                    plot(gxt_shift,rest_scal,'k','linewidth',2);
                    hold on
                    plot(gxt_shift,mod_gxt_shift,'r--','linewidth',2);
                    legend('residual','model');
                    title('shifted/scaled')
                    xx=xlim;

                    subplot(2,2,4)
                    plot(gx,gc,'k','linewidth',2);
                    hold on
                    plot(gx,fullpeak,'r--','linewidth',2);
                    xlim(xx+[-(xx(2)-xx(1)) (xx(2)-xx(1))]+center)
                    ylim(yy);
                end

                %build peak
                ma_global=max(obj.ga.profPF.signal);
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
                pp.peakX_l2r=gx(left:right);
                pp.fullpeak_l2r=fullpeak(left:right);
                pp.paramat=pOpt(:);
                pp.lbound=gx(left);
                pp.rbound=gx(right);
                pp.maxlocation=gx(id);
                pp.height=ma;
                pp.relheight=ma/ma_global;
                pp.modelnumber=1;

                obj.temppeak=peak(pp.peakX_l2r,pp.fullpeak_l2r,pp.modelnumber,pOpt(:),pp);

                %plot temporary peak
                if length(obj.htemp)==2
                    if ishandle(obj.htemp(2))
                        delete(obj.htemp(2));
                    end
                end
                obj.htemp(2)=plot(obj.temppeak,'r',2,obj.mainaxes);
                set(obj.htemp(2),'LineStyle','--')

                disp('add peak - fit: done')


            else
                warning('add peak - fit: htemp(1) is no handle')
            end
        else
            warning('add peak - fit: htemp empty')
        end

    end

    %% ADD PEAK Confirm
    function obj = addPeakOk(obj,hObj,evnt)
        obj.clearHTemp;
        obj.hAdd.btn_addPeakFit.Enabled=false;
        obj.hAdd.btn_addPeakOk.Enabled=false;
        obj.hAdd.btn_addPeakCancel.Enabled=false;
        obj.enableButtons();
        
        if ~isempty(obj.temppeak)
            p=obj.temppeak;
            % get position idx
            cr=find(obj.peakMaxLoc<=p.peakstruct.maxlocation,1,'last');

            % add ga
            obj.ga.peaks=[obj.ga.peaks(1:cr);p; obj.ga.peaks(cr+1:end)];

            %update integrals
            obj.ga.updateIntegrals;

            %update model
            mod=p.getmod(obj.ga.modePF.timegrid);
            obj.ga.modePF.signal=obj.ga.modePF.signal+mod;

            %update model handle
            delete(obj.hModProf);
            obj.hModProf=plot(obj.ga.modePF,'b',1,obj.mainaxes);
            uistack(obj.hModProf,'bottom') %put model plot to the background

            hp=plot(p,'w',2,obj.mainaxes);
            obj.hPeak=[obj.hPeak(1:cr); hp; obj.hPeak(cr+1:end)];
            obj.updateColor;

            %eval fit /overlap
            obj.ga.evalPeaks(max(1,cr-1):min(cr+3,length(obj.ga.peaks)));
            p=obj.ga.peaks(cr+1); % update fit and overlap

            %activate Click CBs in Plot
            obj.actPlotClick

            % update table
            obj.peaktable.set('Visible','off');
            obj.peaktable.addRow({num2str(cr+1,'%i'),num2str(p.peakstruct.maxlocation,'%.3f'),num2str(p.intinfo.abs,'%.3e'), num2str(p.overlap,'%.3f'), num2str(p.fit,'%.3f')},cr+1);
            for i=cr+2:length(obj.ga.peaks)
                obj.peaktable.setCell(i,1,num2str(i,'%i'));
            end
            obj.peaktable.set('Visible','on');

            %update peak locations
            obj.peakMaxLoc=[obj.peakMaxLoc(1:cr); p.peakstruct.maxlocation; obj.peakMaxLoc(cr+1:end)];
            
            obj.peaktable.setSelectedRow(cr+1); % Select correct peak in the table
            obj.addToHistory; % Update history for undo
        else
            warning('add peak - confirm: no peak info available. To be clicked after a peak fit.')
        end
        
    end

    %% ADD PEAK Cancel
    function obj = addPeakCancel(obj,hObj,evnt)
        obj.clearHTemp;
        obj.hAdd.btn_addPeakFit.Enabled=false;
        obj.hAdd.btn_addPeakOk.Enabled=false;
        obj.hAdd.btn_addPeakCancel.Enabled=false;        %obj.hAdd.btn_addPeak.Enabled=true;
        obj.enableButtons();
    end

    %% REFIT EXISTING PEAK
    function obj = ReFitPeak(obj,hObj,evnt)
        if obj.currentSel==0
            errordlg('Select a peak first.')
            return
        end
        % Disable buttons
        isLocked = obj.disableButtons();
        if isLocked
            return
        end
        % Function starts here
        p=obj.ga.peaks(obj.currentSel);
        set(obj.peaktable,'Enable','off');
        obj.deactPlotClick


        %get new fitting
        obj.temppeak=p.refit(obj.ga.profPF.timegrid, obj.ga.profPF.signal, obj.ga.modePF.signal);



        obj.clearHTemp;
        obj.htemp(1)=plot(obj.temppeak,'r',2,obj.mainaxes);
        set(obj.htemp(1),'LineStyle','--')

        disp('add peak - fit: done')


        obj.hMod.btn_PeakReFitOk.Enabled=true;
        obj.hMod.btn_PeakReFitCancel.Enabled=true;
    end

    %% REFIT PEAK OK
    function obj = ReFitOk(obj,hObj,evnt)
        
        obj.clearHTemp;
        obj.hMod.btn_PeakReFitOk.Enabled=false;
        obj.hMod.btn_PeakReFitCancel.Enabled=false;
        obj.enableButtons();

        if ~isempty(obj.temppeak)

            obj.deletePeakRange(obj.currentSel);

            p=obj.temppeak;
            % get position idx
            cr=find(obj.peakMaxLoc<=p.peakstruct.maxlocation,1,'last');
            if isempty(cr) % Catch first peak
                cr=0;
            end
            % add ga
            obj.ga.peaks=[obj.ga.peaks(1:cr);p; obj.ga.peaks(cr+1:end)];

            %update integrals
            obj.ga.updateIntegrals;

            %update model
            mod=p.getmod(obj.ga.modePF.timegrid);
            obj.ga.modePF.signal=obj.ga.modePF.signal+mod;

            %update model handle
            delete(obj.hModProf);
            obj.hModProf=plot(obj.ga.modePF,'b',1,obj.mainaxes);
            uistack(obj.hModProf,'bottom') %put model plot to the background

            hp=plot(p,'w',2,obj.mainaxes);
            obj.hPeak=[obj.hPeak(1:cr); hp; obj.hPeak(cr+1:end)];
            obj.updateColor;

            %eval fit /overlap
            obj.ga.evalPeaks(max(1,cr-1):min(cr+3,length(obj.ga.peaks)));
            p=obj.ga.peaks(cr+1); % update fit, overlap, intinfos

            %activate Click CBs in Plot
            obj.actPlotClick

            % update table
            obj.peaktable.set('Visible','off');

            obj.peaktable.addRow({num2str(cr+1,'%.i'),num2str(p.peakstruct.maxlocation,'%.3f'),num2str(p.intinfo.abs,'%.3e'), num2str(p.overlap,'%.3f'), num2str(p.fit,'%.3f')},cr+1);
            for i=cr+2:length(obj.ga.peaks)
                obj.peaktable.setCell(i,1,num2str(i,'%.i'));
            end
            obj.peaktable.set('Visible','on');

            %update peak locations
            obj.peakMaxLoc=[obj.peakMaxLoc(1:cr); p.peakstruct.maxlocation; obj.peakMaxLoc(cr+1:end)];

            obj.peaktable.setSelectedRow(cr+1); % Select correct peak in the table
            obj.addToHistory; % Update history for undo
        else
            warning('add peak - confirm: no peak info available. To be clicked after a peak fit.')
        end
        obj.actPlotClick
        set(obj.peaktable,'Enable','on');
    end

    %% REFIT PEAK CANCEL
    function obj = ReFitCancel(obj,hObj,evnt)
        obj.clearHTemp;
        obj.hMod.btn_PeakReFitOk.Enabled=false;
        obj.hMod.btn_PeakReFitCancel.Enabled=false;
        obj.enableButtons();
        obj.actPlotClick
        set(obj.peaktable,'Enable','on');
    end

    %% MAINFIG RESIZE FCN
    function sizech_fcn(obj,src,cbd)
        fpos=obj.mainfig.Position;
        obj.peaktable.Units='pixels';
        obj.peaktable.Position=[30 30 330 max(fpos(4)-60,50)];

        obj.mainaxes.Units='pixels';
        obj.mainaxes.Position=[430 60 max(fpos(3)-430-30,50) max(fpos(4)-90,50)];

        set(obj.mainaxes_pos,'Position',[fpos(3)-160 fpos(4)-60 120 18]);
    end

    %% MAINFIG MOUSE MOVE FCN
    function mousemove_fcn(obj,src,cbd)

        loc=get(obj.mainaxes, 'Currentpoint');
        xloc=loc(1,1);
        yloc=loc(1,2);
        xx=get(obj.mainaxes,'XLim');
        yy=get(obj.mainaxes,'YLim');

        % mouse in axes?
        if xloc<xx(1) || xloc>xx(2) || yloc<yy(1) || yloc>yy(2)
            return
        end

        set(obj.mainaxes_pos,'String',[num2str(xloc,'%5.2f') '/' num2str(yloc,'%6.1E')]);
    end

    %% TABLE SELECTION CALLBACK
    function cb_table(obj,hobj,evnt)
        if obj.isSelectCB
            obj.mouseMode([],[],'sel'); % Switch to 'sel', to avoid problems with saving view history
            try
                row=evnt.Indices(1);
                obj.currentSel=row;

                p=obj.ga.peaks(row);
                set(obj.hPeak(row),'linewidth',3)
                rem=1:length(obj.hPeak);
                rem(row)=[];
                set(obj.hPeak(rem),'linewidth',2)

                wi=abs(p.peakstruct.rbound-p.peakstruct.lbound);
                xlim(obj.mainaxes,[p.peakstruct.lbound-0.2*wi p.peakstruct.rbound+0.2*wi]);

                ce_idx=find(obj.ga.profPF.timegrid>=p.peakstruct.maxlocation,1,'first');
                ylim(obj.mainaxes,[0 1.1*obj.ga.profPF.signal(ce_idx)]);


                % split peak info
                obj.hMod.lab_peakidx.Text=num2str(row);
                obj.hMod.lab_submods.Text=num2str(size(p.peakstruct.paramat,2));
            catch
                disp('warning: cell select catch')
            end
        else
            obj.isSelectCB=1;
        end
    end


    %% Create error report
    function obj = errorRep(obj,hObj,evnt)
        save_file = [obj.ga.para.tempf 'ErrorStatus.mat'];
        errordlg(sprintf(['Error-information can be found in the results-folder:\n',...
                [save_file '\n'],...
                'You may contact the developers under:\n',...
                'martina.beese@uni-rostock.de\n',...
                'tomass.andersons@uni-rostock.de\n',...
                'Your request should contain the ErrorLogFile.txt and the ErrorStatus.mat ']));
        obj.ga.errorLog;
    end

end


end

