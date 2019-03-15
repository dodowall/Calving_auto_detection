function [cv_prev,cv_cur,cv_both,mean_cur,mean_both] = cvIm(previousImage,currentImage,segmentationMask,width)
global config;
prev = previousImage;
prev(~segmentationMask) = 0;
cv_prev(1) = coeffVar(double(prev(prev>0)));

cur = currentImage;
cur(~segmentationMask) = 0;
cv_cur(1) = coeffVar(double(cur(cur>0)));

mean_cur(1) = mean(double(cur(cur>0)));

cv_both(1) = coeffVar(abs(double(prev(cur>0)-cur(cur>0))));
mean_both(1) = mean(abs(double(prev(cur>0)-cur(cur>0))));

for i= 1:config.part_num
    prev_i = prev(:,(i-1)*floor(width/config.part_num)+1:i*floor(width/config.part_num));
    cv_prev(i+1) = coeffVar(double(prev_i(prev_i>0)));
    cur_i = cur(:,(i-1)*floor(width/config.part_num)+1:i*floor(width/config.part_num));
    cv_cur(i+1) = coeffVar(double(cur_i(cur_i>0)));
    mean_cur(i+1) = mean(double(cur_i(cur_i>0)));
    
    cv_both(i+1) = coeffVar(abs(double(prev_i(cur_i>0)-cur_i(cur_i>0))));
    mean_both(i+1) = mean(abs(double(prev_i(cur_i>0)-cur_i(cur_i>0))));
end
end