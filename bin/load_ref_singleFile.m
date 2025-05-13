function [dataR] = load_ref_singleFile(ffilename)

% look for existing data
d=dir(ffilename);
if isempty(d)
    error(['filename not found:\n ' ffilename])
end

fol=d.folder;
fname=d.name;

idx=strfind(fname,'.');
if isempty(idx)
    error(['file has no extension: ' fol filesep fname])
end
fname_woE=fname(1:idx(end)-1);

% get all xls and xlsx files in folder that start with the signal file name
% wo extension
fileList = dir(fol);
b=ones(length(fileList),1);

for i=1:length(fileList)
    if fileList(i).isdir
        b(i)=0;
        continue
    end
    
    cfname=fileList(i).name;
    if length(cfname)<4
        b(i)=0;
        continue
    end
    %remove non xls and xlsx
    if ~strcmpi(cfname(end-3:end),'xlsx') && ~strcmpi(cfname(end-2:end),'xls')
        b(i)=0;
        continue
    end
    if ~contains(cfname,'integration','IgnoreCase',1)
        b(i)=0;
        continue
    end
    rep={'rohdaten','Rohdaten','ROHDATEN'};
    fname_woE_rep=replace(fname_woE,rep,'integration');
    if ~startsWith(cfname,fname_woE_rep,'IgnoreCase',1)
        b(i)=0;
        continue
    end
end

if max(b)==0
    error(['determination of reference file failed: ' fol filesep fname])
end
if sum(b)>1
    disp(['multiplie reference files determined: ' fol filesep fname])
    disp(['try detailed comparison of filenames'])
    for i=1:length(fileList)
        if b(i)==0
            continue
        end
        
        cfname=fileList(i).name;
        if length(cfname)<4
            b(i)=0;
            continue
        end

        rep={'rohdaten','Rohdaten','ROHDATEN'};
        fname_woE_rep=replace(fname_woE,rep,'integration');
        if ~strcmp(cfname(length(fname_woE_rep)+1),'.')
            b(i)=0;
            continue
        end
    end
end

if sum(b)>1
    error(['multiplie reference files determined: ' fol filesep fname])
end

fileInfo=fileList(b==1);


xffname=[fileInfo.folder filesep fileInfo.name];
disp('reference data import: found excel reference file')



dataR=xlsread(xffname);

if size(dataR,2)==7
    disp('old reference format detected')
elseif size(dataR,2)==2
    disp('new reference format detected')
else
    error('reference format unknown')
end



























