function [xBoundaries,yBoundaries,calvMask] = change_detection(img1,img2,segmentationMask)
%% change_detection
% Detect changes from a certain area of two consecutive images.
%
%% Syntax
% [xBoundaries,yBoundaries] = change_detection(img1,img2,segmentationMask)

%% Description
% The purpose of this function is to detect textural changes from a certain part of two
% consecutive images. The selected area is represented by the segmentation
% mask. Output of the function are X and Y coordinates of the borders of changed area.

%% Input arguments
% img1 - First image, in grayscale format
% img2 - Second image, captured after the first one. Also in grayscale
% format
% segmentationMask - A binary mask for selecting area to be detected where ones
% represent area of interest. Note that the size of the mask must be equal
% to the size of both two input images

%% Output argument
% xBoundaries - X coordinates of the borders of changed area
% yBoundaries - Y coordinates of the borders of changed area

%% References :
% [1]Y. Zheng, X. Zhang, B. Hou, and G. Liu, Using Combined Difference 
% Image and k-Means Clustering for SAR Image Change Detection, IEEE 
% Geoscience and Remote Sensing Letters, vol. 11, no. 3, pp. 691?695, Mar. 2014.
%
% [2] L. Wang and D.-C. He, Texture classification using texture spectrum, 
% Pattern Recognition, vol. 23, no. 8, pp. 905?910, 1990.
% 
% [3] Rafael C Gonzalez and Richard E Woods. Digital image processing. 
% Dorling Kindersley ; Pearson Prentice Hall, UP, India; 
% [Upper Saddle River, N.J.], 2009.
%

global config;

[height,width]=size(img1);

%% Generate LBP image
nFiltSize=config.LBPFilterSize;
nFiltRadius=config.LBPFilterRadius;
filtR=generateRadialFilterLBP(nFiltSize, nFiltRadius);
textureImg1 = efficientLBP(img1, 'filtR', filtR, 'isRotInv', false, 'isChanWiseRot', false);
textureImg2 = efficientLBP(img2, 'filtR', filtR, 'isRotInv', false, 'isChanWiseRot', false);

textureImg1(~segmentationMask)=0;
textureImg2(~segmentationMask)=0;

textureImg1=double(textureImg1);
textureImg2=double(textureImg2);

%% Calculate two types of difference maps, based on subtraction and log ratio operators
subtractDiff = subtraction(textureImg1,textureImg2);
logDiff = log_ratio(textureImg1,textureImg2);
subtractDiff = uint8(range_normalization(subtractDiff));
logDiff = uint8(range_normalization(logDiff));

%% Apply filters to difference maps
meanFilter = fspecial('average',config.meanFilterSize);
filteredSubtractDiff = imfilter(subtractDiff,meanFilter);
filteredLogDiff = medfilt2(logDiff,config.medianFilterSize);

%% Merge both maps into singel map and threshold change and no-change area
mergedDiff = uint8(config.mergeAlpha * filteredSubtractDiff + (1 - config.mergeAlpha) * filteredLogDiff);
thresholdedMergedDiff = adaptivethreshold(mergedDiff,config.adaptWindow,config.adaptConstant,config.adaptFilterType);

%% Reconstruct polygons from merged map, since changes are detected as clustered points
[yPos,xPos] = find(~thresholdedMergedDiff);
polygons = alphaShape(xPos,yPos,config.alphaShapeRadius,'RegionThreshold',config.alphaShapeRegionThreshold);
xBoundaries = cell(1);
yBoundaries = cell(1);
calvMask = zeros(height,width);

for i=1:numRegions(polygons)
    [~,boundaryCoordinates] = boundaryFacets(polygons,i);
    xBoundaries{i} = boundaryCoordinates(:,1);
    yBoundaries{i} = boundaryCoordinates(:,2);
    calvMask_i = poly2mask(xBoundaries{i},yBoundaries{i},height,width);
    calvMask = calvMask + calvMask_i;
end
end