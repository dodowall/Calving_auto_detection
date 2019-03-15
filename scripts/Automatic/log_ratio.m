function diffImg = log_ratio(img1,img2)
%% log_ratio
% Calculate changes from two images using log ratio operator
%
%% Syntax 
%  diffImg = log_ratio(img1,img2)
%
%% Description
% This function calculate changes from two images using log ratio operator.
% See Equation (2) from [1].

%% Input arguments
% img1 - First image in grayscale format.
% img2 - Second image in grayscale format.

%% Output argument
% diffImg = Difference image
% 

%% Reference
% [1]Y. Zheng, X. Zhang, B. Hou, and G. Liu, ?Using Combined Difference 
% Image and k-Means Clustering for SAR Image Change Detection,? IEEE 
% Geoscience and Remote Sensing Letters, vol. 11, no. 3, pp. 691?695, Mar. 2014.

%% Code
    diffImg = abs(log(img2+1) - log(img1+1));
end