% Test script for runNaturalImageExperiment function
% 
% Notes:
%   1) Run experiment on calibrated 27-in NEC MultiSync PA271 monitor.
%   2) Center observer's eyes horizontal and vertically using chin and forehead rest.
%   3) Set monitor distance to 75 cm.
%   4) Have observer dark adapt in the dark room for 5 min.
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'nPresentations' : (scalar)  Number of presentations per image comparison (default: 5)
%   'controlSignal'  : (string)  Input method for user response (options: 'keyboard', 'gamePad') (default: 'keyboard')
%   'option1Key'     : (string)  For keyboard -> '1', for gamePad either 'GP:UpperLeftTrigger'  or 'GP:X' (default: '1')
%   'option2Key'     : (string)  For keyboard -> '2', For gamePad either 'GP:UpperRightTrigger' or 'GP:A' (default: '2')
%   'giveFeedback'   : (logical) Give feedback if option is on (default: true)
%   'isDemo'         : (logical) Data won't be saved if on (default: false)
%   'subjectName'    : (string)  Name of subject (default: 'testSubject')

runNaturalImageExperiment('experimentName', 'Experiment000',...
                          'nPresentations', 1, ...
                          'controlSignal', 'gamePad', ...
                          'option1Key', 'GP:UpperLeftTrigger', ...
                          'option2Key', 'GP:UpperRightTrigger', ...  
                          'giveFeedback', true, ...
                          'isDemo', false, ...
                          'subjectName', 'AN000');

mglClose;  ListenChar(0);

% NOTE: [q] then [up arrow] to close experiment display