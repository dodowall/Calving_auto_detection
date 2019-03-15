function diff_img = subtraction(img1,img2)
%% subtraction
% Calculate changes from two images using subtraction operator
%
%% Syntax 
%  diffImg = subtraction(img1,img2)
%
%% Description
% This function calculate changes from two images using subtraction operator.
% See Equation (1) from [1].
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
    diff_img = abs(img2 - img1);
end