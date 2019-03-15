function normalized_img = range_normalization(img)
%% range_normalization
% Normalize a matrix value from 0 to 255
%% Syntax
% normalized_img = range_normalization(img)
%% Description
% This function normalize a matrix with arbitraty range into one with range
% from 0 to 255. The normalization process is done using Min-Max
% normalization to form a distribution ranging from 0 to 1 and multiplid by
% 255.

%% Input argument
% img - Input matrix or a grayscale image having values with arbitrary range

%% Output argument
% normalized_img - A matrix or grayscale image with values ranging from 0
% to 255
    min_val = min(img(:));
    max_val = max(img(:));
    range_val = max_val - min_val;
    normalized_img = (img - min_val) / range_val;
    normalized_img = normalized_img * 255;
end