% Test script for analyzeSessionData function
%
% Optional parameters/values:
%   'experimentName' : (string)  Name of experiment folder (default: 'Experiment000')
%   'subjectName'    : (string)  Name of subject (default: 'AN000')
%   'sessionNumber'  : (scalar)  Number of session (default: 1)
%   'plotFigures'    : (logical) Plot figures if option is on (default: true)

analyzeSessionData('experimentName', 'Experiment000', ...
                   'subjectName', 'AN000', ...
                   'sessionNumber', 1, ...
                   'plotFigures', true);