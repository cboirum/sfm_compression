function [lines] = getSegs( I,fillGap,minLength,numPeaks,threshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

BW = edge(I,'canny');
[H,T,R] = hough(BW);

P  = houghpeaks(H,numPeaks,'threshold',ceil(threshold*max(H(:))));
x = T(P(:,2)); y = R(P(:,1));
% Find lines and plot them
lines = houghlines(BW,T,R,P,'FillGap',fillGap,'MinLength',minLength);
max_len = 0;
numLines = length(lines)
for k = 1:numLines
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','blue');

end

