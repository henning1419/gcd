function [gc,gx,ident,flag] = load_gc_csv_singleFile(ffilename)
flag=0;

% determine OS
if isunix
    disp('data import: unix os detected')
end
if ~isunix && ~ismac
    disp('data import: windows os detected')
end


fileList = dir(ffilename);

if isempty(fileList)
    error('data import: check filename. Abort')
else
    disp('data import: found data.')
end

fname=[fileList(1).folder filesep fileList(1).name];
ident=fileList(1).name;

% check file format
isch=0;
iscsv=0;
disp(fname)
if strcmpi(fname(end-1:end),'ch')
    isch=1;
    disp('data import: ch file format detected')
elseif strcmpi(fname(end-2:end),'csv')
    iscsv=1;
    disp('data import: csv file format detected')
end

if ~isch && ~iscsv
    error('data import: check file format. ch and csv not found. Abort')
end



if isch
    temps=ImportAgilentFID(fname);
    if isempty(temps)
        error('data import: could not load ch file')
    end
    if isfield(temps,'Time')
        gx=temps.Time(:);
    else
        error('data import: could not load time grid from ch file')
    end
    if isfield(temps,'Signal')
        gc=temps.Signal(:);
    else
        error('data import: could not load signal from ch file')
    end
    flag=1;
end


if iscsv
    % activate if non utf8 csv files are loaded
    if isunix && ~isdeployed && 0
        conv_sh_path=[cd filesep 'bin' filesep 'convEnc.sh'];
        currp=cd;
        cd(fileList(1).folder)
        system(['sh ' conv_sh_path ]);
        cd(currp)
        fname2=[fname '.utf8'];
    else
        fname2=fname;
    end
    
    ff=fopen(fname2,'r');
    line = fgetl(ff);
    b=0;
    if contains(line,',')
        b=1;
    end
    fclose(ff);
    
    if b
        M=readmatrix(fname2,'filetype','text','delimiter',',');
    else
        M=readmatrix(fname2,'filetype','text');
    end
    % Only on windows without Knime
    % If there are errors, when opening the file, then most likely it is a
    % UTF16-LE file, which is not supported. Then the entries in M are NaN.
    % In this case the file has to be opened with the correct encoding
    % (officially not supported).
    if ~isunix && ~isdeployed && isnan(M(1,1))
        M=readmatrix(fname2,'filetype','text','Encoding','UTF16-LE');
    end
    gx=M(:,1);
    gc=M(:,2);
    
    if isunix && ~isdeployed
        system(['rm "' fname '.utf8"'])
    end
    flag=2;
end















