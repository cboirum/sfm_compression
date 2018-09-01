%% Initialization
reinitialize=0;
if reinitialize
    clear all
    %Color Format:
        colorSpace = 'grayscale'
        %colorSpace = 'rgb'

    adaptorname = 'winvideo';
    dNum = 2;
end
%% loop & display line finding & capture

nFrames = 600;
setFormat =1;
renderVideo = 1;
justOutput = 1; %Only show video that teleoperator will see

%Serial Radio Data Rate Calculations
fRate = 15; %desired fRate for (serial radio)
bitLimiting = 1;% Limit number of points that can be shown (for serial radio)
bitLimit = 115200%9000;% bit/second limit for data output

cornering = 0;
lining = 0;
quadTree = 0;
canny = 0;
cornerHarris = 0;
cornerMinEig = 0;
cornerFAST = 1;
FASTMinContrast = .00125;
FASTMinQuality = .00125;

variableParams = 0;
unionCorners = 0; % 1 = show points that are in either last 2 frames
interCorners = 0; % 1 = show points that are in both last 2 frames
colorSpace = 'RGB';%'grayscale'; 
colorDots = 1;
markerType = 'FilledCircle'; %'Circle'
fps_total = 0;
justVideo = 0;
oldMovie = 0;
regShow = 1;

vidSource = 1;




cornerTest = 0;

if setFormat 
    
    %imageFormat = 'I420_160x120'
    %imageFormat = 'I420_640x480' 

    % imageFormat = 'RGB24_160x120';
     imageFormat = 'RGB24_320x240';
    % imageFormat = 'RGB24_640x480';
    % imageFormat = 'RGB24_1920x1080';
    % imageFormat = 'RGB24_2304x1536';

    %imageFormat = 'MJPG_160x120';
    %imageFormat = 'MJPG_320x240';
    %imageFormat = 'MJPG_640x480';
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
ET=0;
parts = strsplit(imageFormat,'_');
hw = strsplit(parts{2},'x');
W = str2double(hw{1});
H = str2double(hw{2});
overlay = zeros(H,W);
if ~strcmp(colorSpace,'grayscale')
    Red = zeros(H,W);
    Grn = zeros(H,W);
    Blu = zeros(H,W);
    im1c = zeros(H,W);
    bitPerPoint = 44;
else
    bitPerPoint = 28;
end

varSatp = .5;
varSatn = 0.01;
varRng = varSatp-varSatn;
numPnts = 0;

pntLimit = floor(bitLimit/(bitPerPoint*fRate))


for i = 1:nFrames
    
    if variableParams
        fminC = varRng*(i/nFrames)+varSatn;
        fminQ = varRng*(i/nFrames)+varSatn;
        %FASTMinContrast = min(max(fminC,varSatn),varSatp);
        FASTMinQuality = min(max(fminQ,varSatn),varSatp);
    end
    
    if ~fps_total
        tic
    end
    %im1 = rgb2gray(step(vidobj));
    if vidSource == 1
        im1 = step(vidobj);
    elseif vidSource == 2
       
       im1=getdata(vidobj,1,'uint8');
       
    end
    
    if ~strcmp(colorSpace,'grayscale')
       im1c = im1;
       im1 = rgb2gray(im1);
    else
       im1c = im1;
    end
    
    if justOutput
        dispImg = overlay;
    else
        dispImg = im1c;
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
       cornerTest = 1;
    end
    if cornerFAST
       corners = detectFASTFeatures(im1,'MinQuality',FASTMinQuality,'MinContrast',FASTMinContrast);
       cornerTest = 1;
    end
    if cornerMinEig
       cornersAll = detectMinEigenFeatures(im1);
       cornerTest = 1;
    end
    if cornerTest
       if bitLimiting
            corners = corners.selectStrongest(pntLimit);
       end
       npnts = [corners.Location,3+zeros(size(corners.Location,1),1)];
       if unionCorners
           pnts = union(lastpnts,npnts,'rows');
           lastpnts = npnts;
       elseif interCorners
           pnts = intersect(lastpnts,npnts,'rows');
           lastpnts = npnts;
       else
           pnts = npnts;
       end
       numPnts = numPnts +length(pnts);

       if colorDots
           xy = pnts(:,1:2);
           x = xy(:,2);
           y = xy(:,1);
           indx = sub2ind(size(im1),x',y');
           if strcmp(colorSpace,'grayscale')
               % get intensity of image at dot location 
               colors = im1(indx);
               colors = [colors' colors' colors'];
           else
               %get color of image at each dot location
               Red = im1c(:,:,1);
               Grn = im1c(:,:,2);
               Blu = im1c(:,:,3);
               r = Red(indx);
               g = Grn(indx);
               b = Blu(indx);
               colors = [r' g' b'];
           end
           
           dispImg = insertShape(dispImg,markerType,pnts,'Color',colors,'Opacity',1);
       else
           dispImg = insertShape(dispImg,markerType,pnts);
       end
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
    
    fps_ave_s = sprintf('average fps: %.1f',(i/ET));
    if renderVideo
        
        if regShow
            try
            fps_s = sprintf('        fps: %.1f',fps);
            dispImg = insertText(dispImg,[10,10],fps_s);
            dispImg = insertText(dispImg,[10,30],fps_ave_s);
            catch
            end

        end
        if ~oldMovie
            hvpc.step(dispImg);
        end
    
    
    end
    if ~fps_total
        t = toc;
        fps = 1/t;
        
        ET = ET+t;%;(fps_ave+fps)/i;
        
    end
    %hold off
end

imageFormat;
averagePoints = numPnts/i;
FASTMinQuality;
FASTMinContrast;
fps_ave_s;
aveBitrate = averagePoints*(bitPerPoint*fRate);
fprintf('%s \t ave bit/sec: %.1f  \t ave points: %.1f \t min qual: %.3f \t min cont: %.3f \t Format: %s Rendered: %i \n',fps_ave_s,aveBitrate,averagePoints,FASTMinQuality,FASTMinContrast,imageFormat,renderVideo)
if vidSource ==2
    stop(vidobj);
end
if fps_total
    t = toc;
    fps = nFrames/t;
end
% %% clear out memory
% 
% release(vidobj);
% clear vidobj;