function bicomb_wrapper

addpath functions % make files in the /functions/ folder accessible

obsInitials  = input('\nEnter subject initials:', 's');
infoFilePath = sprintf('rawdata/%s_ObsInfo.mat', obsInitials);

dFile = sprintf('diaries/%s-%s.txt',obsInitials,date);
diary(dFile)

if exist(infoFilePath, 'file')
    load(infoFilePath, 'OB')
else
    yes = input('Make new subject? (y/n)', 's');
    if strcmp(yes,'y') % yes
        OB = defobserver(obsInitials);
    else               % no
        error('No subject for those initials')
    end
end

EX = defexperiment;

bicomb_core(EX,OB)

rmpath functions % remove /Functions/ folder from the current path

return

function EX = defexperiment

load('monitor/amonitor.mat','MON'); % monitor file, currently a dummy
EX = MON;
EX.whichScreen   = 0;
EX.patchDuration = input('\nPatching duration for this test session (in minutes)?:');

EX.isBefore = 0;
befAftStr = '';
while (~strcmp(befAftStr,'A')) && (~strcmp(befAftStr,'B'))
    befAftStr = input('\nIs this before/after patching? (B/A):', 's');
end

EX.isBefore = strcmp(befAftStr,'B');
EX.baseContrast = 0.5; % 50% base contrast
EX.nTrials      = 8; 

if EX.isBefore
    EX.timeSincePatch = 0;
    EX.contrastRatio  = [0.5, 1, 2]; % half, equal, double
else
    EX.timeSincePatch = input('\nNumber of minutes since patch was removed?:');
    EX.contrastRatio  = input('\nContrast ratio to test?:');
end

return

function OB = defobserver(obsInitials)

    OB.initials = obsInitials;
    OB.fullName = input('\nEnter subject name:', 's');
    OB.patchEye = '';
    while (~strcmp(OB.patchEye,'L')) && (~strcmp(OB.patchEye,'R'))
        OB.patchEye = input('\nWhich eye will be patched? (L/R):', 's');
    end
    
    OB.notes = input('\nEnter any notes:', 's');
    OB.dataPath = sprintf('rawdata/%s/', OB.initials);
    mkdir(OB.dataPath)
    
    OB.infoFilePath = sprintf('rawdata/%s_ObsInfo.mat', OB.initials);
    OB.eyeStr = {'L','R'};
                          
    save(OB.infoFilePath, 'OB')

return