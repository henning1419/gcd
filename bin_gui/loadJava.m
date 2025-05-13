function varargout = loadJava()


%% Add java class paths

% Disable warnings
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame')
warning('off','MATLAB:Java:DuplicateClass');
warning('off','MATLAB:javaclasspath:jarAlreadySpecified')

filename = mfilename('fullpath');
filename = regexprep(filename, 'gQSPsim/app/(\.\./app)+', 'gQSPsim/app');
filename = strrep(filename , '/../app', '');
RootPath = fileparts(filename);
% Tree controls and XLWRITE
if isdeployed
    %     Paths = {
    %         'UIExtrasTree.jar';
    %         'UIExtrasTree.jar'
    %         };
else
    Paths = {
        fullfile(RootPath,'+uix','+resource','UIExtrasTable.jar')
        fullfile(RootPath,'+uix','+resource','UIExtrasTree.jar')
        };
    
    % Add paths
    javaaddpath(Paths);
end



end
