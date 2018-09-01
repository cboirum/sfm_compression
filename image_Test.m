%% CV Project: Simgle Image Test version


try
    oops = firstTime
catch
    firstTime  = 1
end

if firstTime
%vlFeat required for fast gaussian convolution:
run('C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\vlfeat-0.9.18\toolbox\vl_setup')
end

% reinitialize=0;
% if reinitialize
%     clear all
%     %Color Format:
%         colorSpace = 'grayscale'
%         %colorSpace = 'rgb'
% 
%     adaptorname = 'winvideo';
%     dNum = 1;
% end
%% loop & display line finding & capture

fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\fig1.png'
fPath = 'C:\Users\BOIRUM\Documents\Curt\Research Ideas\Video Polygon Compression\Lenna.png'
%fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\Test images\nasa_mars_rover_video.jpg';

%fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\Test images\Picture1.png';
%fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\fig1.png'
%fPath = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\Test images\newYork1.png'
outDir = 'C:\Users\BOIRUM\Documents\Curt\CMU\Fall 2013\Computer Vision\Project\Test images\';
testOut = 1;
writeOutput = 0;
dNum = 2;
nFrames = 1;
setFormat =0;
renderVideo = 1;
justOutput = 1; %Only show video that teleoperator will see
edgeSeed = 1;
xSeeds = 4;
ROI_AVG_Coloring = 1;

%Serial Radio Data Rate Calculations
fRate = 30; %desired fRate for (serial radio)
fRateEnforce = 0; %Enforce display rate
bitLimiting = 1;% Limit number of points that can be shown (for serial radio)
%bitLimit = 115200%9000;% bit/second limit for data output
bitLimit = 200000


%Items to Display
dotImage = 1; %Display image of feature points/verices 
upscale = 1; %upscale display resolution
upscale2 = 1; %scale uniformily
scaleFactor = .5; 
% uW = 640;
% uH = 480;
uH = 600;
uW = 800;

cornerFAST = 1;
FASTMinContrast = .000125;
FASTMinQuality = .0125;

%Extras
cornering = 0;
lining = 0;
quadTree = 0;
canny = 0;
cornerHarris = 0;
cornerMinEig = 0;


showTriangles = 1;


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
dotSize = 1; %size of vertex markers

cDotSize = 3; %size of center markers
triangulate = 1; %Use delaunay triangulation to join corners
plotTriangle = 0; %plot triangulation externally
centerType = 'centroid'; %type of triangle center to use for fill color lookup
cColorType = 'squareMean'; %use square of radius r to calc average color of triangle

%gradient = 1; %calculate gradient

gSmooth = 1; %Gaussian smoothing before processing
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

im1 = imread(fPath);
hw = size(im1(:,:,1));
W = hw(2);
H = hw(1);
overlay = zeros(H,W);
vertexImage = zeros(H,W);
if upscale
    rW = uH/H;
    rH = uW/W;
end
if upscale2
    rW = scaleFactor;
    rH = scaleFactor;
    uH = floor(H*rH);
    uW = floor(W*rW);
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
    plot(boundsY,boundsX,'o')
    %[boundsX';boundsY'];
else
    boundsX = [1 1 H H]';
    boundsY = [1 W 1 W]';
end
tic
for i = 1:nFrames
    
    if variableParams
        fminC = varRng*(i/nFrames)+varSatn;
        fminQ = varRng*(i/nFrames)+varSatn;
        %FASTMinContrast = min(max(fminC,varSatn),varSatp);
        FASTMinQuality = min(max(fminQ,varSatn),varSatp);
    end
    
    
    
    if ~strcmp(colorSpace,'grayscale')
       im1c = im1;
       im1 = rgb2gray(im1);
    else
       im1c = im1;
    end
    
    if gSmooth
       im1 = vl_imsmooth(im2double(im1),gSigma);
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
    dispImgRIO_AVG = dispImg;
    
    if cornering
        tic
        corners = corner(im1);
        tCornerFunc = toc
        cornerTest = 1
        %plot(corners(:,1), corners(:,2), 'r*');
        %regShow = 1;
    end
    if cornerHarris
        tic
       corners = detectHarrisFeatures(im1,'MinQuality', 0.1);
       tHarris = toc
       cornerTest = 1;
    end
    if cornerFAST
        tic
       corners = detectFASTFeatures(im1,'MinQuality',FASTMinQuality,'MinContrast',FASTMinContrast);
       tFAST = toc
       cornerTest = 1;
    end
    if cornerMinEig
        tic
       cornersAll = detectMinEigenFeatures(im1);
       tMinEig = toc
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
           x1 = x;
           y = xy(:,1);
           y1 = y;
           if triangulate
               
               tri = delaunay(double(x),double(y));
               if plotTriangle
                    triplot(tri,y,x);
               end
               
               polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];
               polyPntsUS = polyPnts;
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
           
           indx = sub2ind(size(im1),x',y');
%            if strcmp(cColorType,'squareMean')
%                color = imAveColSqr(im1c, cx, cy, r, overlay);
%            end
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
           if ROI_AVG_Coloring
               %Use average color of each triangle/polygon instead of just
               %color of central pixel. (polyPntsUS = unscaled points)

               ROI_AVG_Colors = getRoiAveColors(im1c,polyPntsUS);
           end
           %Use gausian of triangle area (perhaps using polygon circle as
           %radius for gaussian? Variable Gaussian radias could be a
           %problem. 
           
           if showDots
                dispImg = insertShape(dispImg,markerType,pnts,'Color',colors,'Opacity',1);
           end
           if dotImage
               %clear vertexImage
               vertexImage(:) = uint8(0);
               %vertexImage = insertShape(vertexImage,markerType,polyPnts,'Color',colors,'Opacity',1);
               vertexImage = insertMarker(uint8(vertexImage),[y1,x1],'size',1,'Color',colors(1:length(x1),:));
               %vertexImage = insertMarker(im1c,[y,x],'size',1,'Color',colors);
               vertexImage2 = insertMarker(im1c,[y1,x1],'o','color','g','size',1);
               vertexImage3 = insertMarker(uint8(vertexImage),[pcY,pcX],'size',1,'Color',pcColors);
               imshow(vertexImage2)
           end
           if triangulate
                if showTriangles
                    
                    dispImg = insertShape(dispImg,'FilledPolygon',polyPnts,'Color',pcColors,'Opacity',1);
                    if ROI_AVG_Coloring
                        dispImgRIO_AVG = insertShape(dispImg,'FilledPolygon',polyPnts,'Color',ROI_AVG_Colors,'Opacity',1);
                    end
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
        tic
        lines= getSegs(im1,fillGap,minLength,numPeaks,threshold);
        tHough = toc
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
    dispImg = uint8(dispImg);
    if justVideo
        dispImg = im1;
    end
    if canny
        regShow = 0;
        BW = edge(im1,'canny');
        %imshow(BW)
        dispImg = BW;
    end
tTotal = toc
subplot(2,2,1), subimage(dispImg);

subplot(2,2,2), subimage(im1c);

subplot(2,2,3), subimage(uint8(dispImgRIO_AVG))

end
% 
% imageFormat;
% averagePoints = numPnts/i;
% FASTMinQuality;
% FASTMinContrast;
% fps_ave_s;
% aveBitrate = averagePoints*(bitPerPoint*fRate);
% fprintf('%s \t ave bit/sec: %.1f  \t ave points: %.1f \t min qual: %.3f \t min cont: %.3f \t Format: %s Rendered: %i \n',fps_ave_s,aveBitrate,averagePoints,FASTMinQuality,FASTMinContrast,imageFormat,renderVideo)
numTriangles
% if vidSource ==2
%     stop(vidobj);
% end
% if fps_total
%     t = toc;
%     fps = nFrames/t;
% end
% %% clear out memory
% 
% release(vidobj);
% clear vidobj;
firsTime = 0;
if testOut
    numVerts = length(x);
    numTri = length(tri);

    %Minimum typcasting
    triMin = uint16(tri);
    xMin = uint16(x);
    yMin = uint16(y);
    pcColorsMin = uint8(pcColors);

    [path,inName,ext] = fileparts(fPath);
    
    %fname = sprintf('testIms/testIMG_VectCompSml2_fps%d_bLim%d.mat',[fRate,bitLimit]);
    fname = [outDir, inName, sprintf('_outPut_fps%d_bLim%d',[fRate,bitLimit]),'.mat'];
    save(fname,'xMin','yMin','triMin','pcColorsMin')
    fname2 = [outDir, inName, sprintf('_outPut_fps%d_bLim%d',[fRate,bitLimit]),'.bmp'];
    imwrite(dispImg,fname2,'bmp');
    fname3 = [outDir, inName, sprintf('_outPut_fps%d_bLim%d',[fRate,bitLimit]),'.jpg'];
    imwrite(dispImg,fname3,'jpg');
end
% %% Save output in different formats
% if writeOutput
%     fName = 'testIms/testIMG_';
%     formats = imformats
%     labels = {};
%     sizes = [];
%     for i = 1:length(formats)
%         good = 1;
%         fmt = formats(i).ext{1}
%         try
%             name = [fName,fmt,'.',fmt];
%         imwrite(im1c,name,fmt);
%         catch
%             fprintf('bad format: %s\n',fmt)
%             bad = 1;
%         end      
%     end
% end
% 
% %% binary encoding that requires rebuilding polygon points using:
% %polyPnts = [y(tri(:,1)),x(tri(:,1)),y(tri(:,2)),x(tri(:,2)),y(tri(:,3)),x(tri(:,3)),y(tri(:,1)),x(tri(:,1))];
% 
% if writeOutput
%     numVerts = length(x)
%     numTri = length(tri)
% 
%     %Minimum typcasting
%     triMin = uint16(tri);
%     xMin = uint16(x);
%     yMin = uint16(y);
%     pcColorsMin = uint8(pcColors);
% 
%     fname = sprintf('testIms/testIMG_VectCompSml2_fps%d_bLim%d.mat',[fRate,bitLimit]);
%     save(fname,'xMin','yMin','triMin','pcColorsMin')
%     fname2 = sprintf('testIms/testIMG_VectCompSml2_fps%d_bLim%d.bmp',[fRate,bitLimit]);
%     imwrite(dispImg,fname2,'bmp');
% end