function J = bilateralFCol(I, windowSize, sigmaSpatial,sigmaRange)
I = double(I)/255;
ILab = rgb2lab(I);
m = size(I,1);
n = size(I,2);
J = zeros(m,n);
hN = (windowSize-1)/2;
%Initializing the gaussian
gaussSpatial = zeros(windowSize, windowSize);
for kDiff = -hN:hN
    for lDiff = -hN:hN
        gaussSpatial(hN+kDiff+1,hN+lDiff+1) = exp(-(((kDiff)^2 + (lDiff)^2)/(2*sigmaSpatial^2)));
    end
end
% LAB VALUES RANGE FROM [0,100]
% at first the sigmaRange value was left unscaled but the image produced
% was not affected by the filter. It was later discovered by arbritarily
% setting sigmaRange to large values that the sigma was too small.
% The sigma is therefore scaled by the largest LAB value
sigmaRange = 100*sigmaRange;


% The bilateral filter formulation is based on the discrete formulation of
% the bilateral filter dicussed in the report.
for i = 1:m
    for j = 1:n
        %Boundary Handling - The window and gaussian are bounded, values
        %which lie outside the image bounds but within the window are 
        %ignored. (padding with zeros produces black edges on the image)
    
        %Defining Bounds
        % Only conisdering values that lie withing the image and window
        xLowerBound = max(i-hN,1); xUpperBound = min(i+hN,m);
        yLowerBound = max(j-hN,1); yUpperBound = min(j+hN,n);
        windowL = ILab(xLowerBound:xUpperBound, yLowerBound:yUpperBound,1);
        windowA = ILab(xLowerBound:xUpperBound, yLowerBound:yUpperBound,2);
        windowB = ILab(xLowerBound:xUpperBound, yLowerBound:yUpperBound,3);
        
        %Gaussian Boundary Handling - 
        % The spatial gauss function has to be resized to match the 
        % windows bounded dimensions 
        gaussXBound = (xLowerBound:xUpperBound)+hN+1-i;
        gaussYBound = (yLowerBound:yUpperBound)+hN+1-j;
        
        %Computing the photometric gaussian
        squareDiffRangeL = (ILab(i,j,1)-windowL).^2;
        exponentRange   = squareDiffRangeL/(2*sigmaRange^2);
        gaussRangeL      = exp(-exponentRange);
        
        squareDiffRangeA = (ILab(i,j,2)-windowA).^2;
        exponentRange   = squareDiffRangeA/(2*sigmaRange^2);
        gaussRangeA      = exp(-exponentRange);
        
        squareDiffRangeB = (ILab(i,j,3)-windowB).^2;
        exponentRange   = squareDiffRangeB/(2*sigmaRange^2);
        gaussRangeB      = exp(-exponentRange);
        
        %Computing the weight
        WL = gaussRangeL.*gaussSpatial(gaussXBound,gaussYBound);
        WA = gaussRangeA.*gaussSpatial(gaussXBound,gaussYBound);
        WB = gaussRangeB.*gaussSpatial(gaussXBound,gaussYBound);
        
        %Returning the normalzied function 
        % The bilateral filter is the product of the 
        J(i,j,1) = sum(WL(:).*windowL(:))/sum(WL(:));
        J(i,j,2) = sum(WA(:).*windowA(:))/sum(WA(:));
        J(i,j,3) = sum(WB(:).*windowB(:))/sum(WB(:));
    end
end
J = lab2rgb(J);
