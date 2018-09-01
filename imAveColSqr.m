function [ color ] = imAveColSqr( I, cx, cy, r, overlay )
%imAveColSqr returns average color or intensity for square of half-length r
%and centered at cx, cy

Isize = size(I);
w = Isize(1);
h = Isize(2);

L = cx - r;
R = cx + r;
T = cy - r;
B = cy + r;

roi = overlay;
roi(T:B,L:R) = 1;
color = regionprops(roi, I, 'MeanIntensity');

% if size(size(I))>2
%     %color image
% %     Red = mean(mean(I(L:R,T:B,1)));
% %     Grn = mean(mean(I(L:R,T:B,2)));
% %     Blu = mean(mean(I(L:R,T:B,3)));
% %     color = [Red,Grn,Blu];
%     
% else
%     %grayscale image
% %     color = mean(mean(I(L:R,T:B)));
% end

end

