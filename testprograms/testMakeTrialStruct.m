%%testMakeTrialStruct test script for makeTrialStruct function

makeTrialStruct('directoryName','RandomTargetShapeFixedIlluminantFixedBkGnd',...
    'LMSstructName', 'LMSStruct',...
    'outputFileName', 'exampleTrial',...
    'nBlocks', 2,...
    'stdYIndex', 5, ...
    'cmpYIndex', (1:10));