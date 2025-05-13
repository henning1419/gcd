function [dataRc] = load_ref(folname)
flag=0;

% look for existing data

fileList = dir([folname filesep '*.xlsx']);


if isempty(fileList)
    disp('data import: no mat or csv data found in data folder. Abort')
    return
else
    disp('reference data import: found xlsx data')
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
dataRc=cell(length(fileList),1);



for i=1:length(fileList)
    dataRc{i}=xlsread([fileList(i).folder filesep fileList(i).name]);
    
end
flag=1;































