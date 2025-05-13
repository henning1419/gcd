function [Dc,Tc,flag] = load_gc_csv(folname)
flag=0;

% look for existing data
if isfile([folname filesep 'data.mat'])
    disp('data import: data.mat exists in data folder - using this one')
    load([folname filesep 'data'],'Dc','Tc');
    flag=1;
else
    fileList = dir([folname filesep '*.csv']);
    
    if isempty(fileList)
        fileList = dir([folname filesep '*.CSV']);
    else
        disp('data import: found csv data')
    end
    
    if isempty(fileList)
        disp('data import: no mat or csv data found in data folder. Abort')
        return
    else
        disp('data import: found CSV data')
    end
    
    
    
    % sort files by minute after last _ in file name
    strnames=fileList;
    numnames=[];
    for i=1:length(strnames)
        pos = regexp(strnames(i).name, '_');
        strnames(i).name=strnames(i).name(pos(end)+1:end-4);
        numnames=[numnames str2num(strnames(i).name)];
    end
    
    [~,idx]=sort(numnames);
    fileList=fileList(idx);
    
    
    % load data from csv
    Tc={};
    Dc={};
    
    % handle UTF-16LE format (works only under unix)
    if isunix
        conv_sh_path=[cd filesep 'bin' filesep 'convEnc.sh'];
        currp=cd;
        cd(fileList(1).folder)
        system(['sh ' conv_sh_path ]);
        cd(currp)
        
        for i=1:length(fileList)
            M=readmatrix([fileList(i).folder filesep fileList(i).name '.utf8'],'filetype','text');
            Tc{end+1}=M(:,1);
            Dc{end+1}=M(:,2);
            system(['rm "' fileList(i).folder filesep fileList(i).name '.utf8"'])
        end
        
    else
        for i=1:length(fileList)
            M=readmatrix([fileList(i).folder filesep fileList(i).name],'filetype','text');
            Tc{end+1}=M(:,1);
            Dc{end+1}=M(:,2);
        end
        
    end
    
    
    save([fileList(i).folder filesep 'data.mat'], 'Tc','Dc');
    flag=1;
end












