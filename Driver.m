I = imread('lena.tif');
%I = rgb2gray(I);
windowSize = 7;
sigmaSpatial = 15;
sigmaRange = 0.1;
Ifilter = bilateralFGray(I,windowSize,sigmaSpatial,sigmaRange);


figure
imshow(I)
figure
imshow(Ifilter)