function [registeredImage,tForm,prev] = registration4(descriptorsReference,descriptorsLocationReference,imageToRegister,tForm_previous)
%% registration
% Register an input image with a reference image. 
%% Syntax
% registeredImage = registration(descriptorsReference,descriptorsLocationReference,imageToRegister)
%
%% Description
% This function is intended to geometrically align the input image to the
% reference image. The reference image is represented by its descriptors
% and descriptors' location.

%% Input arguments
% descriptorReference - Descriptors of the reference image
% descriptorsLocationReference - Location of the descriptors of the
% reference image
% imageToRegister - Grayscale image to be registered.

%% Output argument
% registeredImage - Aligned version of the input image
%% References
% [1] Herbert Bay, Andreas Ess, Tinne Tuytelaars, and Luc Van Gool. 
% Speeded-up robust features (surf). Comput. Vis. Image Underst., 
% 110(3):346?359, June 2008.
% [2] Z. Xiong and Y. Zhang, ?A critical review of image registration 
% methods,? International Journal of Image and Data Fusion, vol. 1, 
% no. 2, pp. 137?158, Jun. 2010.
% [3] P. H. S. Torr and A. Zisserman, ?MLESAC: A New Robust Estimator with
% Application to Estimating Image Geometry,? Computer Vision and Image 
% Understanding, vol. 78, no. 1, pp. 138?156, Apr. 2000.


%% Extract SURF feature from  the input image
surfInput = detectSURFFeatures(imageToRegister);
%% Extract descriptors and their corresponding location from SURF feature of
%input image
[descriptorsInput,descriptorsLocationInput] = extractFeatures(imageToRegister,surfInput);
%% Find matched descriptors from input and reference images
indexPairs = matchFeatures(descriptorsReference,descriptorsInput);
if length(indexPairs)>=10 && length(descriptorsLocationInput.Location)>=200
    %% Select location of matched features
    matchedPointsReference = descriptorsLocationReference(indexPairs(:,1));
    matchedPointsInput = descriptorsLocationInput(indexPairs(:,2));
    %% Estimate geometric transform and noise removal
    [tForm,~,~] = estimateGeometricTransform(matchedPointsInput,matchedPointsReference,'affine');
    prev = 0;
else
    tForm = tForm_previous;
    prev = 1;
end
outputView = imref2d(size(imageToRegister));
registeredImage = imwarp(imageToRegister,tForm,'OutputView',outputView);
end