%% loop & display line finding & capture
colorSpace = 'rgb'
adaptorname = 'winvideo';

dNum = 2;
nFrames = 200;
fRate = 30; %desired fRate for (serial radio)
setFormat =1;
renderVideo = 1;
justOutput = 0; %Only show video that teleoperator will see

%Serial Radio Data Rate Calculations

fRateEnforce = 0; %Enforce display rate

vidSource = 1;

if setFormat 
    
    %imageFormat = 'I420_160x120'
    %imageFormat = 'I420_640x480' 

    % imageFormat = 'RGB24_160x120';
    % imageFormat = 'RGB24_320x240';
    % imageFormat = 'RGB24_640x480';
    % imageFormat = 'RGB24_1920x1080';
    % imageFormat = 'RGB24_2304x1536';

    %imageFormat = 'MJPG_160x120';
    %imageFormat = 'MJPG_320x240';
    imageFormat = 'MJPG_640x480';
    %imageFormat = 'MJPG_1920x1080';
    try
        release(vidobj)
    catch
    end
    try
        stop(vidobj)
    catch
    end
    if vidSource == 1
        vidobj = imaq.VideoDevice(adaptorname, dNum);
        imaqInfo = imaqhwinfo;
        hwInfo = imaqhwinfo(adaptorname);
        hwInfo.DeviceInfo(dNum);
        device1 = hwInfo.DeviceInfo(dNum);
        set(vidobj, 'VideoFormat', imageFormat);
        set(vidobj, 'ReturnedColorSpace', colorSpace);
    elseif vidSource == 2
        vidobj = videoinput(adaptorname, dNum, imageFormat);
        src = getselectedsource(vidobj);
        vidobj.FramesPerTrigger =1;
        vidobj.TriggerRepeat = Inf;
        vidobj.ReturnedColorspace = colorSpace; 
        if fRateEnforce
            src.FrameRate = fRate;
        else
            src.FrameRate = '30';
        end
    end
end
%%
%Try better performance video player for speedup
hvpc = vision.VideoPlayer;

if vidSource ==2
    fprintf('Starting device')
    start(vidobj)
end


for i = 1:nFrames
    if vidSource == 1
        im1 = step(vidobj);
    elseif vidSource == 2
        im1=getdata(vidobj,1,'uint8');
    end
    if renderVideo
        hvpc.step(im1);
    end
    pause(1/fRate);
    outPath = sprintf('testpic%d.jpg',i)
    imwrite(im1,outPath);
end