function NaturalImageThresholdsLocalHook
% NaturalImageThresholdsLocalHook
%
% Configure things for working on the NaturalImageThresholds project.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localHookFolder) and delete "Template" from the filename,
% this will get run when you execute tbUseProject('NaturalImageThresholds') 
% to set up for this project.
%
% You will need to edit your local copy for your computer's configuration.

%% Say hello
theProject = 'NaturalImageThresholds';
fprintf('Running %s local hook\n',theProject);

%% Remove old preferences
if (ispref(theProject))
    rmpref(theProject);
end

%% Set path for Dropbox
sysInfo = GetComputerInfo();
switch (sysInfo.localHostName)
    case 'eagleray'
        % DHB's office desktop
        baseDir = fullfile(filesep,'Volumes','Users1','Dropbox (Aguirre-Brainard Lab)');
        calFile = 'NaturalImageThresholdsCal_Philly';
    case 'Amys-iMac-Pro'
        % Amy's desktop
        baseDir = fullfile('/Users','amyni','Dropbox (Aguirre-Brainard Lab)');
        calFile = 'NaturalImageThresholdsCal_Amy';
    otherwise
        % Some unspecified machine, try user specific customization
        switch(sysInfo.userShortName)
            % Could put user specific things in, but at the moment generic
            % is good enough.
            otherwise
                baseDir = fullfile('/Users',sysInfo.userShortName,'Dropbox (Aguirre-Brainard Lab)');
        end
        calFile = 'NaturalImageThresholdsCal_Philly';
end

%% Set base directory for project materials
BaseDir = fullfile(baseDir,'CNST_materials',theProject);
setpref(theProject,'BaseDir',BaseDir);

%% Calibration folder/file
setpref(theProject,'CalDataFolder',fullfile(BaseDir,'CalibrationData'));
setpref('BrainardLabToolbox','CalDataFolder',fullfile(BaseDir,'CalibrationData'));
setpref(theProject,'CalDataFile',calFile);

%% Device string for computer that runs this project's calibration
setpref('BrainardLabToolbox','PR650DevicePortString','/dev/cu.KeySerial1');

%% End