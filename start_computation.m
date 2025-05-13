% Specify the folder with .csv or .ch files (without the last / or \)
% This folder can only contain these files, other files will throw an error
% In this folder a new folder for results will be created
tempf=[cd filesep 'example_data'];

% Specify a folder to copy the results to
savefol=[cd filesep 'lastCalc' filesep];

% prep and run
folderInfo = dir(tempf);% read elements of the folder
mkdir(tempf,'results') % create results folder
for i=1:length(folderInfo)
    if folderInfo(i).isdir==0
        main([tempf filesep folderInfo(i).name],'calc');
        copyfile([savefol 'gca_res.mat'],[tempf filesep 'results' filesep folderInfo(i).name '.mat']) % copy and rename results file to folder
        copyfile([savefol 'peak_res.csv'],[tempf filesep 'results' filesep folderInfo(i).name 'peak_res.csv']) % copy and rename results file to folder
    end
end
