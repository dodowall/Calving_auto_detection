% load config file. Config contains parameters for some methods
load config;
%define config as global variable
global config;

%get a list of all image files on selected path
imageFiles = dir(fullfile(imagePath,strcat('*.',config.imageFormat)));
startImageNumber = str2double(regexp(startImageFilename,'(\d*)','match'));
endImageNumber = str2double(regexp(endImageFilename,'(\d*)','match'));

%filter the images to contain a set of selected sequence, from start number
%to end number
fileNamesToFilter = strcat(config.imagePrefix,'_',arrayfun(@num2str,(startImageNumber:endImageNumber)','UniformOutput', false),'.',config.imageFormat);
imageFiles = imageFiles(ismember({imageFiles.name},fileNamesToFilter));

imageFileDate = {imageFiles.date}';
imageFileName = {imageFiles.name}';

% load satellite data reference
load_Href = load ('H.mat');
Href = load_Href.H;

% if ~exist('currentIm.mat')
     lastit = 0;
% else
%     load('calv_auto.mat');
% 	load('calv_part_auto.mat');
% 	load('calv_coord_auto.mat');
%     load('calv_cv_image_auto.mat');
%     load('calv_cv_both_auto.mat');
%     load('initial.mat');
%     load('currentIm.mat');
% end

for i=lastit+1:min(lastit+1+nit,length(imageFileName))
    if i==1

        [~,imName,~] =  fileparts(imageFileName{i});
        calvName_i{i}=strjoin({imName,num2str(i),datestr(datenum(imageFileDate(i)),'ddmmyyyy_HHMM')},'_');
        cv_image{i} =  [];
        cv_both{i} = [];        
        currentImage_i{i} = [];
        currentMask_i{i} = [];
            
        load('initial.mat')

    else
        previousImage = currentImage;
        previousMask = currentMask;
        
        currentImage0 = imread(fullfile(imagePath,imageFileName{i}));
        currentImage0 = rgb2gray(currentImage0);
        
        %register the current image using the first image as reference
        currentImage0 = registration(descriptorsReference,descriptorsLocationReference,currentImage0,in);
        currentImage = imcrop(currentImage0,[169 1189 3911 443]);
        currentImage_i{i} = currentImage;
        
        % segment the current image based on the initial mask selected on
        % the first image.
        currentMask = activecontour(currentImage,initialMask,config.activeContourIteration);
        currentMask_i{i} = currentMask;
        
        %combine the previous segmentation mask and current one to form
        %combined segmentation mask
        segmentationMask = previousMask & currentMask;
                
        % Calculate the coefficient of variation
        [cv_prev_i,cv_cur_i,cv_both_i] = cvIm(previousImage,currentImage,segmentationMask,width);
        
        cv_image{i} =  mean([cv_cur_i' cv_prev_i']');
        cv_both{i} = cv_both_i;
        [~,imName,~] =  fileparts(imageFileName{i});
        calvName_i{i}=strjoin({imName,num2str(i),datestr(datenum(imageFileDate(i)),'ddmmyyyy_HHMM')},'_');
        disp(strcat(calvName_i{i},' done')); 

    end
end

%saving the data

calving_currentImage_auto = cell2struct(currentImage_i,calvName_i,2);
calving_currentMask_auto = cell2struct(currentMask_i,calvName_i,2);
calving_cv_image_auto = cell2struct(cv_image,calvName_i,2);
calving_cv_both_auto = cell2struct(cv_both,calvName_i,2);
save('calv_current_auto.mat','calving_currentImage_auto','calving_currentMask_auto');
save('calv_cv_auto.mat','calving_cv_image_auto','calving_cv_both_auto');
