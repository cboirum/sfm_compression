%% cam info
% set resetCamera = 1 if running for first time since Matlab boot

% resetCamera = 1
% 
% 
% vidSource = 1; %use faster video streaming source
% 
% if resetCamera && vidSource == 1
%     clear all
%     vidSource = 1;
%     %dNum = 3; % default webcam
%     dNum = 2; %c920
%     imaqInfo = imaqhwinfo;
%     hwInfo = imaqhwinfo('winvideo');
%     hwInfo.DeviceInfo(dNum)
%     device1 = hwInfo.DeviceInfo(dNum)
%     defaultFormat = device1.DefaultFormat
%     supportedFormats = device1.SupportedFormats
% 
%     %Cam Format:
%     
%     %imageFormat = 'I420_160x120'
%     imageFormat = 'I420_640x480' 
% 
%     % imageFormat = 'RGB24_160x120'
%     % imageFormat = 'RGB24_640x480';
%     % imageFormat = 'RGB24_1920x1080';
%     % imageFormat = 'RGB24_2304x1536';
% 
%     %imageFormat = 'MJPG_160x120';
%     %imageFormat = 'MJPG_320x240';
%     %imageFormat = 'MJPG_640x480';
%     %imageFormat = 'MJPG_1920x1080';
% 
%     ind=find(ismember(supportedFormats,imageFormat));
%     if isempty(ind)
%         fprintf('Unsupported image format, using default instead:\n')
%         imageFormat = defaultFormat
%     end
% 
%     fprintf('Image format: %s\n',imageFormat)
% 
%     %Color Format:
%     colorSpace = 'grayscale'
% %     colorSpace = 'rgb'
% 
%     %%% Configure Camera
% 
%     adaptorname = 'winvideo';
% 
%     deviceid = dNum;
% 
%     clear vidobj
% 
%     vidobj = imaq.VideoDevice(adaptorname, dNum);
%     release(vidobj)
%     set(vidobj, 'VideoFormat', imageFormat);
%     set(vidobj, 'ReturnedColorSpace', colorSpace);
%     
% end
% % %% GUI display of video feed
% % clear vidobj
% % vidobj = videoinput('winvideo',dNum, imageFormat);
% % figure('Toolbar','none',...
% %        'Menubar', 'none',...
% %        'NumberTitle','Off',...
% %        'Name','Camera Feed Preview');
% % 
% % 
% % set(vidobj, 'ReturnedColorSpace', colorSpace);   
% % vidRes = get(vidobj, 'VideoResolution')
% % nBands = get(vidobj, 'NumberOfBands');
% % hImage = image( zeros(vidRes(2), vidRes(1), nBands) );
% %    
% % preview(vidobj, hImage);
% % 
% % 
% % 
% 
% % %% Time capture of multiple images
% % tic
% % frames = {}
% % clear im1
% % for i = 1:1
% %     frames{i} = step(vidobj);
% %     
% % end
% % toc
% % im1 = frames{i};
% % %% preview last image
% % figure
% % imshow(im1)
% % 
% % %% get line segments in image
% % imshow(im1)
% % hold on
% % getSegs(im1)
% % hold off
% %% Alternate Video Sourcing Method
% if vidSource == 2
%     vidobj = videoinput('winvideo', 2, 'MJPG_160x120'); %select input device
% 
%     hvpc = vision.VideoPlayer;   %create video player object
% 
%     src = getselectedsource(vidobj);
%     vidobj.FramesPerTrigger =1;
%     vidobj.TriggerRepeat = Inf;
%     %vidobj.ReturnedColorspace = 'rgb';
%     vidobj.ReturnedColorspace = 'grayscale';
%     src.FrameRate = '30';
%     try
%         start(vidobj)
%     catch
%     end
% 
%     %start main loop for image acquisition
%     for t=1:500
%       imgO=getdata(vidobj,1,'uint8');    %get image from camera
%       hvpc.step(imgO);    %see current image in player
%     end
% end
clear all
%Color Format:
    colorSpace = 'grayscale'
    %colorSpace = 'rgb'

adaptorname = 'winvideo';
dNum = 2;
%% loop & display line finding & capture
regShow = 1;
nFrames = 180
cornering = 0;
lining = 0;
quadTree = 0;
canny = 0;
cornerHarris = 0;
cornerMinEig = 0;
cornerFAST = 1;
FASTMinContrast = .1;
FASTMinQuality = .1;
fps_total = 0;
justVideo = 0;
oldMovie = 0;
setFormat =1;
vidSource = 1;
cornerTest = 0;

if setFormat 
    
    %imageFormat = 'I420_160x120'
    %imageFormat = 'I420_640x480' 

    % imageFormat = 'RGB24_160x120'
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
        vidobj = videoinput(adaptorname, dNum, imageFormat)
        src = getselectedsource(vidobj);
        vidobj.FramesPerTrigger =1;
        vidobj.TriggerRepeat = Inf;
        vidobj.ReturnedColorspace = colorSpace; 
        src.FrameRate = '30';
    end
end

%Try better performance video player for speedup
if ~oldMovie
    hvpc = vision.VideoPlayer;
end
if vidSource ==2
    fprintf('Starting device')
    start(vidobj)
end
if fps_total
    tic
end
fps_ave=0;

for i = 1:nFrames
    if ~fps_total
        tic
    end
    %im1 = rgb2gray(step(vidobj));
    if vidSource == 1
        im1 = step(vidobj);
    elseif vidSource == 2
       
       im1=getdata(vidobj,1,'uint8');
       dispImg = im1;
    end
    %imshow(im1)
    %hold on
    if cornering
        corners = corner(im1);
        cornerTest = 1
        %plot(corners(:,1), corners(:,2), 'r*');
        %regShow = 1;
    end
    if cornerHarris
       corners = detectHarrisFeatures(im1,'MinQuality', 0.1);
       valid_corners = corners;
       %%[features, valid_corners] = extractFeatures(im1, corners);
%        RGB = insertShape(RGB,'FilledPolygon',Pts_poly.',...
%                               'Color',[0 1 1],'Opacity',0.2); 
%        RGB = insertShape(RGB,'Line',TwoLanes1',...
%             'Color',{'yellow','magenta'});
       corners = valid_corners;
       cornerTest = 1;
       %plot(valid_corners);
    end
    if cornerFAST
       corners = detectFASTFeatures(im1,'MinQuality',FASTMinQuality,'MinContrast',FASTMinContrast);
       %corners = detectFASTFeatures(im1);
       %corners = corners.selectStrongest(100);
       cornerTest = 1;
    end
    if cornerMinEig
       cornersAll = detectMinEigenFeatures(im1);
       cornerTest = 1;
    end
    if cornerTest
       pnts = [corners.Location,3+zeros(size(corners.Location,1),1)];
       im1 = insertShape(im1,'Circle',pnts);
    end
    if lining
        regShow = 1;
        fillGap = 1;
        minLength = 2;
        numPeaks = 5;
        threshold = 0.2;
        lines= getSegs(im1,fillGap,minLength,numPeaks,threshold);
        numLines = length(lines)
    end
    if quadTree
       regShow = 0;
       %im1 = imread('liftingbody.png');
       [H,W] = size(im1);
       im1Pad = padarray(im1,[(512-H)/2,(512-W)/2]);
       %im1Crop = im1(1:128, 1:128);
       S = qtdecomp(im1Pad,.9);
       blocks = repmat(uint8(0),size(S));
       for dim = [512 256 128 64 32 16 8 4 2 1];
           numblocks = length(find(S==dim));
           if numblocks >0
               values = repmat(uint8(1),[dim dim numblocks]);
               values(2:dim,2:dim,:) = 0;
               blocks = qtsetblk(blocks,S,dim,values);
           end
       end
       imshow(blocks,[])
    end
    if justVideo
        dispImg = im1;
    end
    if canny
        regShow = 0;
        BW = edge(im1,'canny');
        %imshow(BW)
        dispImg = BW;
    end
    if regShow
        try
        fps_s = sprintf('        fps: %.1f',fps);
        fps_ave_s = sprintf('average fps: %.1f',((i-1)/fps_ave));
        im1 = insertText(im1,[10,10],fps_s);
        im1 = insertText(im1,[10,30],fps_ave_s);
        catch
        end
        dispImg = im1;
    end
    if ~oldMovie
        hvpc.step(dispImg);
    else
        imshow(dispImg);
        ffff = 100
        M(i) = getframe;
    end
    
    
    
    if ~fps_total
        t = toc;
        fps = 1/t;
        
        fps_ave = fps_ave+t;%;(fps_ave+fps)/i;
        
    end
    %hold off
end
if vidSource ==2
    stop(vidobj)
end
if fps_total
    t = toc
    fps = nFrames/t
end
% %% clear out memory
% 
% release(vidobj);
% clear vidobj;