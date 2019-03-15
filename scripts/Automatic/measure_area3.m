function [calvingSize,calvingSizePart,ncalv] = measure_area3(xBoundaries,yBoundaries,xQ,yQ,actualPixelSize)
%% measure_area
% Measure the area in real size from sets of pixels area represented by
% their borders. 

%% Syntax
% calvingSize = measure_area(xBoundaries,yBoundaries,xQ,yQ,hReference)

%% Description
% This function 'converts' pixel size of the area of changes from pixel to
% real size (meter) based on data derived from a satellite.

%% Input arguments
% xBoundaries - X coordinates of the borders of changed area
% yBoundaries - Y coordinates of the borders of changed area
% xQ - 
% yQ - 

%% Output argument
% calvingSize - Array contains area for each detected change in square meter.
global config;
    
% Width of image
[height,width] = size(xQ);
calvArea = zeros(height,width);
calvingSizePart = zeros(config.part_num+1,1);
ncalv = zeros(config.part_num+1,1);

    %% In case no change is detected, return empty value
    if isempty(xBoundaries{1})
        calvingSize={[]};
    else
        %% Calculate the area in real size for each detected change
        calvingSize = cell(1);
        ncalv(1) = length(xBoundaries);
        for i = 1:length(xBoundaries)
            in = poly2mask(xBoundaries{i},yBoundaries{i},height,width);%inpolygon(xQ,yQ,xBoundaries{i},yBoundaries{i});
            calvArea(in) = actualPixelSize(in);
            calvingSize{i} = sum(calvArea(in));
            for j=1:config.part_num
                in_j = in(:,(j-1)*floor(width/config.part_num)+1:j*floor(width/config.part_num));
                if sum(in_j(:))>config.alphaShapeRegionThreshold
                    ncalv(j+1) = ncalv(j+1)+1;
                end
            end
                    
        end
        calvingSizePart(1) = sum(calvArea(:));
        for j=1:config.part_num
            calvArea_j = calvArea(:,(j-1)*floor(width/config.part_num)+1:j*floor(width/config.part_num));
            calvingSizePart(j+1) = sum(calvArea_j(:));  
        end
    end
end