function J = bilateralFGray(I, windowSize, sigmaSpatial,sigmaRange)
I = double(I)/255;
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
        window = I(xLowerBound:xUpperBound, yLowerBound:yUpperBound);
        
        % Gaussian Boundary Handling - 
        % The spatial gauss function has to be resized to match the 
        % windows bounded dimensions 
        gaussXBound = (xLowerBound:xUpperBound)+hN+1-i;
        gaussYBound = (yLowerBound:yUpperBound)+hN+1-j;
        
        %Computing the photometric gaussian
        squareDiffRange = (I(i,j)-window).^2;
        exponentRange   = squareDiffRange/(2*sigmaRange^2);
        gaussRange      = exp(-exponentRange);
        
        %Computing the weight
        W = gaussRange.*gaussSpatial(gaussXBound,gaussYBound);
        
        %Returning the normalzied function 
        J(i,j) = sum(W(:).*window(:))/sum(W(:));
    end
end
