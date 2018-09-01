%% cam info
clear all
dNum = 3; % default webcam
dNum = 2; %c920
imaqInfo = imaqhwinfo;
hwInfo = imaqhwinfo('winvideo');
hwInfo.DeviceInfo(dNum)
device1 = hwInfo.DeviceInfo(dNum)
defaultFormat = device1.DefaultFormat
supportedFormats = device1.SupportedFormats

%Cam Format:


% imageFormat = 'RGB24_640x480';
% imageFormat = 'RGB24_1920x1080';
% imageFormat = 'RGB24_2304x1536';

imageFormat = 'MJPG_160x120';
%imageFormat = 'MJPG_320x240';
%imageFormat = 'MJPG_640x480';
%imageFormat = 'MJPG_1920x1080';

ind=find(ismember(supportedFormats,imageFormat));
if isempty(ind)
    fprintf('Unsupported image format, using default instead:\n')
    imageFormat = defaultFormat
end

fprintf('Image format: %s\n',imageFormat)

%Color Format:
colorSpace = 'grayscale'

%%% Configure Camera

adaptorname = 'winvideo';

deviceid = dNum;

clear vidobj

vidobj = imaq.VideoDevice(adaptorname, dNum);

set(vidobj, 'ReturnedColorSpace', colorSpace);
set(vidobj, 'VideoFormat', imageFormat);

%% GUI display of video feed
clear vid
vid = videoinput('winvideo',dNum, imageFormat);
figure('Toolbar','none',...
       'Menubar', 'none',...
       'NumberTitle','Off',...
       'Name','Camera Feed Preview');


set(vid, 'ReturnedColorSpace', colorSpace);   
vidRes = get(vid, 'VideoResolution')
nBands = get(vid, 'NumberOfBands');
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );
   
preview(vid, hImage);




%% Time capture of multiple images
tic
frames = {}
clear im1
for i = 1:1
    frames{i} = step(vidobj);
    
end
toc
im1 = frames{i};
%% preview last image
figure
imshow(im1)

%% get line segments in image
imshow(im1)
hold on
getSegs(im1)
hold off
%% loop & display line finding & capture
tic
nFrames = 90
cornering = 0;
lining = 0;
quadTree = 1;
for i = 1:nFrames
    im1 = rgb2gray(step(vidobj));
    %im1 = step(vidobj);
    %imshow(im1)
    hold on
    if cornering
        corners = corner(im1);
        plot(corners(:,1), corners(:,2), 'r*');
    end
    if lining
        fillGap = 1;
        minLength = 2;
        numPeaks = 5;
        threshold = 0.2;
        lines= getSegs(im1,fillGap,minLength,numPeaks,threshold);
        numLines = length(lines)
    end
    if quadTree
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
    end
    imshow(blocks,[])
    M(i) = getframe;
    hold off
end
t = toc
fps = nFrames/t
%% clear out memory

release(vidobj);
clear vidobj;