function bicomb_core(EX,OB)

ST  = defstim(EX);
rep = 1;
switch EX.isBefore
    case 0
        fStem = sprintf('%s_patchdur%0.0f_before', ...
                        OB.initials, EX.patchDuration);
    case 1
        fStem = sprintf('%s_patchdur%0.0f_after_%0.0f', ...
                        OB.initials, EX.patchDuration, EX.timeSincePatch);
end

resultsFile = sprintf('%s/%s_r%0.0f.mat',OB.dataPath,fStem,rep);

while exist(resultsFile,'file')
    rep = rep + 1;
    resultsFile = sprintf('%s/%s_r%0.0f.mat',OB.dataPath,fStem,rep);
end

RS.nT = zeros(2,length(EX.contrastRatio));
RS.linePosition = nan(EX.nTrials,2,length(EX.contrastRatio));

try

    % PSYCHTOOLBOX INITIALISATION
    [windowPtr, destRect] = initpsychtoolbox(EX);
    
    while any(isnan(RS.linePosition(:)))
        iPolarity = randi(2);
        while ~any(isnan(reshape(squeeze(RS.linePosition(:,iPolarity,:)).',1,[])))
            iPolarity = randi(2);
        end

        if length(EX.contrastRatio) > 1
            iContrastRatio = randi(length(EX.contrastRatio));
            while ~any(isnan(reshape(RS.linePosition(:,iPolarity,iContrastRatio).',1,[])))
                iContrastRatio = randi(length(EX.contrastRatio));
            end
        else
            iContrastRatio = 1;
        end

        iTrial = RS.nT(iPolarity,iContrastRatio) + 1;

        switch iPolarity
            case 1
                eyePhases = [-22.5, 22.5];
            case 2
                eyePhases = [22.5, -22.5];
        end

        for iEye = 2:-1:1
            ST.phaseDeg = eyePhases(iEye);
            imGrat(:,:,iEye) = getgratstim(ST);
        end

        ST.lineY = randi(round((ST.gratHDeg/2)*10))/10 - ST.gratHDeg/4;

        isConfirmed = 0;
        
        i = 0;

        while ~isConfirmed
            
            showstimulus(windowPtr, imGrat, EX, ST)

            isKeyDown = 0;
            lineShift = 0.1;

            while ~isKeyDown
                [~, ~, keyCode] = PsychHID('KbCheck');
                if keyCode(KbName('ESCAPE'))
                    error('Escape key pressed, experiment aborted.')
                elseif keyCode(KbName('UpArrow'))
                    ST.lineY  = ST.lineY - lineShift;
                    isKeyDown = 1;
                elseif keyCode(KbName('DownArrow'))
                    ST.lineY  = ST.lineY + lineShift;
                    isKeyDown = 1;
                elseif keyCode(KbName('SPACE'))
                    isConfirmed = 1;
                    isKeyDown   = 1;
                end
            end
            
            waitSecs(0.01)
            
        end

        RS.linePosition(iTrial,iPolarity,iContrastRatio) = ST.lineY;
        RS.nT(iPolarity,iContrastRatio) = RS.nT(iPolarity,iContrastRatio) + 1;

        ST.overlay   = getfusionlock(ST); % get new fusion pattern before new trial

    end

catch ME
    
    Screen('CloseAll')
    save('lasterror.mat')
    
end

save(resultsFile,'RS','EX','OB')

screen('CloseAll')

return

function showstimulus(windowPtr, imGrat, EX, ST)

    imLines  = zeros(ST.imSz);
    lineYPix = round(ST.lineY * ST.pixPerDeg);
    xOn = ST.imSz/4; xOff = xOn+round(ST.imSz/8);
    yOn = ST.imSz/2-1+lineYPix; yOff = yOn+2;
    imLines(yOn:yOff,xOn:xOff) = -1;
    imLines = imLines + fliplr(imLines);
    
    r1 = [0 0 EX.scrHeight EX.scrHeight];
    destRect = CenterRectOnPoint(r1, EX.scrWidth*0.5, EX.scrHeight*0.5);
    
    for i = 0:1

        imToShow = (1+imGrat(:,:,i+1) + ST.overlay + imLines)./2;

        Screen('SelectStereoDrawBuffer', windowPtr, i);
        Screen('FillRect', windowPtr, 0.5);
        stimFrame = Screen('MakeTexture', windowPtr,  imToShow, [], [], 2);
        Screen('DrawTexture', windowPtr, stimFrame, [], destRect);

    end
    Screen('Flip', windowPtr);

return

function [windowPtr, destRect] = initpsychtoolbox(EX)

    LoadPsychHID;
    KbName('UnifyKeyNames');
    
    PsychImaging('PrepareConfiguration')
    PsychImaging('AddTask', 'General', 'SideBySideCompressedStereo')
    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
    PsychImaging('AddTask', 'General', 'EnablePseudoGrayOutput');
    
    [windowPtr,~] = PsychImaging('OpenWindow', EX.whichScreen, [], [], [], 102);
    r1 = [0 0 EX.scrHeight EX.scrHeight];
    destRect = CenterRectOnPoint(r1, EX.scrWidth*0.5, EX.scrHeight*0.5);

return

function ST = defstim(EX)

    ST.imSz         = EX.scrHeight; % (pixels)
    ST.pixPerDeg    = EX.pixPerDeg; 
    ST.gratWDeg     = 6;   % (deg)
    ST.gratHDeg     = 6;
    ST.spatSigmaPix = 2;   % (pix)
    ST.fusGrain  = ST.pixPerDeg/2;
    ST.fusBrd    = ST.pixPerDeg;
    ST.thetaDeg  = 90;
    ST.spatFreq  = 0.3;
    
    if ST.imSz/ST.fusGrain ~= round(ST.imSz/ST.fusGrain)
        error('ST.imSz/ST.fusGrain ~= round(ST.imSz/ST.fusGrain)')
    end

    ST.overlayC  = 0.9;
    ST.overlay   = getfusionlock(ST);

return

function imGrat = getgratstim(ST)

    thetaRad = ST.thetaDeg * (pi/180);
    phaseRad = ST.phaseDeg * (pi/180);
    x = linspace(-(ST.imSz/ST.pixPerDeg)*pi, ...
                  (ST.imSz/ST.pixPerDeg)*pi,ST.imSz);

    [xx,yy] = meshgrid(x);
    u = xx .* cos(thetaRad) - yy .* sin(thetaRad);

    imSin = sin(u.*ST.spatFreq+phaseRad);

    imEnvW = ones(ST.imSz);
    imEnvH = ones(ST.imSz);
    
    x = linspace(-(ST.imSz/(2*ST.pixPerDeg)), ...
                  (ST.imSz/(2*ST.pixPerDeg)),ST.imSz);
    [xx,yy] = meshgrid(x);
    
    imEnvW(xx<(-ST.gratWDeg/2)) = 0;
    imEnvW(xx>(ST.gratWDeg/2+1/ST.pixPerDeg)) = 0;
    imEnvH(yy<(-ST.gratHDeg/2)) = 0;
    imEnvH(yy>(ST.gratHDeg/2+1/ST.pixPerDeg)) = 0;
    
    imGauss = getgauss(64,ST.spatSigmaPix);
    
    imEnv = conv2(imEnvH .* imEnvW, imGauss, 'same');

    imGrat = imSin .* imEnv;

return

function im = getgauss(imSize,sigmaPix)

    x = (1:imSize)-imSize/2-0.5;
    [xx,yy] = meshgrid(x,x);

    im = normal(xx, sigmaPix) .* normal(yy, sigmaPix);

return

function y = normal(x,sigma)

    y = (1/(sigma*sqrt(2*pi))) .* exp(-x.^2/(2*sigma^2));

return

function imLock = getfusionlock(ST)

    noiseSamples = randi(2,[ST.imSz/ST.fusGrain,ST.imSz/ST.fusGrain])*2-3;
    imNoise      = imresize(noiseSamples,ST.fusGrain,'nearest');
    imMask       = zeros(ST.imSz);
    imMask((3*ST.fusBrd+1):end-3*ST.fusBrd,(3*ST.fusBrd+1):end-3*ST.fusBrd) = 1;
    imMask((4*ST.fusBrd+1):end-4*ST.fusBrd,(4*ST.fusBrd+1):end-4*ST.fusBrd) = 0;
    imLock = imNoise .* imMask .* 1;

return