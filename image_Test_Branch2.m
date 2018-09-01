%% Initialization
%required for fast gaussian convolution:
run('C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\vlfeat-0.9.18\toolbox\vl_setup')

reinitialize=1;
if reinitialize
    clear all
    %Color Format:
        colorSpace = 'grayscale'
        %colorSpace = 'rgb'

    adaptorname = 'winvideo';
    dNum = 1;
end
%% loop & display line finding & capture

fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\Test images\Picture1.png';
im2 = imread(fPath);

fileInput = 0; %Use the above input file instead of video feed

H2 = size(im2,1);
W2 = size(im2,2);

dNum = 2;
nFrames = 100;
setFormat =1;
renderVideo = 1;
justOutput = 1; %Only show video that teleoperator will see

%Serial Radio Data Rate Calculations
fRate = 30; %desired fRate for (serial radio)
fRateEnforce = 0; %Enforce display rate
bitLimiting = 1;% Limit number of points that can be shown (for serial radio)
%bitLimit = 115200%9000;% bit/second limit for data output
bitLimit = 200000

upscale = 1; %upscale display resolution
% uW = 640;
% uH = 480;
uH = 600;
uW = 800;


cornering = 0;
lining = 0;
quadTree = 0;
canny = 0;
cornerHarris = 0;
cornerMinEig = 0;
cornerFAST = 1;
FASTMinContrast = .000125;
FASTMinQuality = .0125;

showTriangles = 1;
edgeSeed = 0;
showCenterDots = 0;
if justOutput
    separateOriginal = 0; %Separate video of raw incoming stream
else
    separateOriginal = 1;
end
separateVerts = 0; %Video of verticies
separateCenters = 0; %Video of Centers

showFPS = 0; %Overlay fps
showDots = 0; %Render dots
dotSize = 3; %size of vertex markers

cDotSize = 3; %size of center markers
triangulate = 1; %Use delaunay triangulation to join corners
plotTriangle = 0; %plot triangulation externally
centerType = 'centroid'; %type of triangle center to use for fill color lookup
%gradient = 1; %calculate gradient

gSmooth = 0; %Gaussian smoothing before processing
gSigma = .5; %Gaussian kernal sigma

variableParams = 0;
unionCorners = 0; % 1 = show points that are in either last 2 frames
interCorners = 0; % 1 = show points that are in both last 2 frames
colorSpace = 'RGB';%'grayscale'; %

colorDots = 1;
markerType = 'FilledCircle'; %'Circle'
fps_total = 0;
justVideo = 0;
oldMovie = 0;
regShow = .5;

vidSource = 1;




cornerTest = 0;

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

%Try better performance video player for speedup
if ~oldMovie
    hvpc = vision.VideoPlayer;
    hvpc2 = vision.VideoPlayer;
end
if vidSource ==2
    fprintf('Starting device')
    start(vidobj)
end
if fps_total
    tic
end
ET=0;
if fileInput
    W = W2;
    H = H2;
else
    parts = strsplit(imageFormat,'_');
    hw = strsplit(parts{2},'x');
    W = str2double(hw{1});
    H = str2double(hw{2});
end
overlay = zeros(H,W);
if upscale
    rH = uH/H;
    rW = uW/W;
end
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

imSizeEstBITS = bitPerPoint*pntLimit

numTriangles = 0;
if edgeSeed
    %Definition of image border seed points for equal distribution to reduce 
    %edge distortion of triangles
    xSeeds = 20;
    seedDist = floor(W/xSeeds);
    xTopBot = [1:seedDist:W,W];
    yLeftRight = [1:seedDist:H,H];
    xCount = length(yLeftRight);
    yCount = length(xTopBot);
    
    xLeft = 1+zeros(1,xCount);
    xRight = W + zeros(1,xCount);
    yTop = 1 + zeros(1,yCount);
    yBot = H + zeros(1,yCount);

    boundsY = [xTopBot(1:end-1) xLeft(2:end) xRight xTopBot(2:end-1)]';
    boundsX = [yTop(1:end-1) yLeftRight(2:end) yLeftRight yBot(2:end-1)]';
    %plot(boundsY,boundsX,'o')
    %[boundsX';boundsY'];
else
    boundsX = [1 1 H H]';
    boundsY = [1 W 1 W]';
end

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
    if fileInput
        im1 = im2;
    else
        if vidSource == 1
            im1 = step(vidobj);
        elseif vidSource == 2

           im1=getdata(vidobj,1,'uint8');

        end
    end
    
    if gSmooth
       im1 = vl_imsmooth(im1,gSigma);
    end
    
    if ~strcmp(colorSpace,'grayscale')
       im1c = im1;
       im1 = rgb2gray(im1);
    else
       im1c = im1;
    end
    
    if justOutput
        
        if upscale
            dispImg = zeros(uH,uW);
        else
            dispImg = overlay;
        end
    else
        dispImg = im1c;
    end

    
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
       npnts = [corners.Location,dotSize+zeros(size(corners.Location,1),1)];
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
           
           pnts = [pnts;[boundsY boundsX dotSize + zeros(length(boundsX),1)]];
%            x2 = [x;[1 1 H H]'];
%            y2 = [y;[1 W 1 W]'];
           xy = pnts(:,1:2);
           x = xy(:,2);
           y = xy(:,1);
           if triangulate
               
               tri = delaunay(double(x),double(y));
               if plotTriangle
                    triplot(tri,y,x);
               end
               
               polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];
               if strcmp(centerType,'centroid')
                   yIndx = [1 3 5];
                   xIndx = [2 4 6];
                   polyCntrs = [ceil(mean(polyPnts(:,yIndx),2)),ceil(mean(polyPnts(:,xIndx),2))];
                   if upscale
                       xIndx = [2 4 6 8];
                       yIndx = [1 3 5 7];
                       polyPnts(:,xIndx) = floor(polyPnts(:,xIndx)*rW);
                       polyPnts(:,yIndx) = floor(polyPnts(:,yIndx)*rH);
                   end
                   
               elseif strcmp(centerType,'incenter')
               end
%                x2 =[x;polyCntrs(:,2)];
%                y2 =[y;polyCntrs(:,1)];
               pcX = polyCntrs(:,2);
               pcY = polyCntrs(:,1);
               numTriangles = length(pcY);
               x2 =[boundsX;pcX];
               y2 =[boundsY;pcY];
               pcCount = length(pcX);
               pnts = [pnts;[y2 x2 dotSize+zeros(length(x2),1)]];
               x =[x;x2];
               y =[y;y2];
               
           end
           
           %cntrIndx = sub2ind(size(im1),pcX',pcY');
           
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
           pcColors = colors(end-pcCount+1:end,:);
           if showDots
                dispImg = insertShape(dispImg,markerType,pnts,'Color',colors,'Opacity',1);
           end
           
           if triangulate
                if showTriangles
                    
                    dispImg = insertShape(dispImg,'FilledPolygon',polyPnts,'Color',pcColors,'Opacity',1);
                end
                if showCenterDots
                    
                    dispImg = insertShape(dispImg,'FilledCircle',[pcY,pcX,dotSize + zeros(length(pcX),1)],'Color',pcColors,'Opacity',1);
                end
           end
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
        
        if showFPS
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
        if separateOriginal
            hvpc2.step(im1c);
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
numTriangles
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

%% Save output in different formats

fName = 'testIms/testIMG_';
formats = imformats
labels = {};
sizes = [];
for i = 1:length(formats)
    good = 1;
    fmt = formats(i).ext{1}
    try
        name = [fName,fmt,'.',fmt];
    imwrite(im1c,name,fmt);
    catch
        fprintf('bad format: %s\n',fmt)
        bad = 1;
    end      
end

%% Encode vectorized image

% %2d polygon points in M x 8 array
% %polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];
% polyPnts;
% %8bit colors in M x 3 array (RGB) or M x 1 (grayscale)
% pcColors;
% im1cVect = struct();
% im1cVect.polyList = polyPnts;
% im1cVect.polyColors = pcColors;
% save('testIms/testIMG_VectComp.mat','im1cVect')

% %% 2d polygon points in M x 7 array
% %polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];
% polyPnts2 = polyPnts(:,1:7);
% %8bit colors in M x 3 array (RGB) or M x 1 (grayscale)
% pcColors;
% im1cVect = struct();
% im1cVect.polyList = polyPnts2;
% im1cVect.polyColors = pcColors;
% save('testIms/testIMG_VectComp2.mat','im1cVect')

%% binary encoding that requires rebuilding polygon points using:
%polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];

numVerts = length(x)
numTri = length(tri)

%Minimum typcasting
triMin = uint16(tri);
xMin = uint16(x);
yMin = uint16(y);
pcColorsMin = uint8(pcColors);

fname = sprintf('testIms/testIMG_VectCompSml2_fps%d_bLim%d.mat',[fRate,bitLimit]);
save(fname,'xMin','yMin','triMin','pcColorsMin')
fname2 = sprintf('testIms/testIMG_VectCompSml2_fps%d_bLim%d.bmp',[fRate,bitLimit]);
imwrite(dispImg,fname2,'bmp');