function calving_quantification(imagePath,startImageFilename,endImageFilename,nit,yearFolder,resultFolder)
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
load_Href = load ([yearFolder,'H.mat']);
Href = load_Href.H;

if ~exist([resultFolder,'currentIm.mat'])
    lastit = 0;
else
    load([resultFolder,'calv_auto.mat']);
	load([resultFolder,'calv_part_auto.mat']);
	load([resultFolder,'calv_coord_auto.mat']);
    load([resultFolder,'calv_cv_image_auto.mat']);
    load([resultFolder,'calv_cv_both_auto.mat']);
    load([resultFolder,'initial.mat']);
    load([resultFolder,'currentIm.mat']);
end

for i=lastit+1:min(lastit+1+nit,length(imageFileName))
    if i==1
        calving_i{i} = [];
        coord_i{i} = [];
        [~,imName,~] =  fileparts(imageFileName{i});
        calvName_i{i}=strjoin({imName,num2str(i),datestr(datenum(imageFileDate(i)),'ddmmyyyy_HHMM')},'_');
        calving_part_i{i} = [];
        cv_image{i} =  [];
        cv_both{i} = [];
        currentImage_i{i} = [];
        currentMask_i{i} = [];
        
        if ~exist([resultFolder,'initial.mat'])
            %read image and convert to grayscale
            currentImage0 = imread(fullfile(imagePath,imageFileName{i}));
            currentImage0 = rgb2gray(currentImage0);

            [height0,width0]=size(currentImage0);

            %extract SURF features, descriptors and descriptors' location of
            %the first image in a sequence
            surfReference = detectSURFFeatures(currentImage0);
            [descriptorsReference,descriptorsLocationReference] = extractFeatures(currentImage0,surfReference);

             %show the image and let user select ocean
            figure;
            imshow(currentImage0);
            [xocean,yocean]=getline;
            close;    

            initialoceanMask = poly2mask(xocean,yocean,height0,width0);
            oceanMask0 = activecontour(currentImage0,initialoceanMask,config.activeContourIteration);
            [x_ocean,y_ocean] = mask2poly(oceanMask0,0,0);        

            figure;
            imshow(currentImage0)
            hold on
            plot(x_ocean,y_ocean)
            uiwait(gcf);

            %Crop image, the cropped size is matched the benchmark image size
            currentImage = imcrop(currentImage0,config.cropRect);
            oceanMask = imcrop(oceanMask0,config.cropRect);
            [height,width]=size(currentImage);

            hReference  = meshgrid(Href,zeros(height,1));
            [xQ,yQ]=meshgrid(1:width,1:height);

            %show the image and let user select an initial segmentation mask
            figure;
            imshow(currentImage);
            [x,y]=getline;
            close;

            %segment the image usign Chan-vese algorithm based on the initial
            %mask
            initialMask = poly2mask(x,y,height,width);
            currentMask = activecontour(currentImage,initialMask,config.activeContourIteration);
            % Remove points situated in the ocean
            currentMask(oceanMask == 1) = 0;

            % Points situated at the front and in the ocean
            [x_poly,y_poly] = mask2poly(currentMask,config.cropRect(1), config.cropRect(2));        

            figure;
            imshow(currentImage0)
            hold on
            plot(x_poly,y_poly)
            plot(x_ocean,y_ocean)
            uiwait(gcf);

            currentMask0 = poly2mask(x_poly,y_poly,height0,width0);
            mask0 = currentMask0 | oceanMask0;
            [x0,y0] = mask2poly(mask0,0,0);          
            in = inpolygon(descriptorsLocationReference.Location(:,1),descriptorsLocationReference.Location(:,2),x0,y0);

            save([resultFolder,'initial.mat'],'initialMask','currentImage','currentMask',...
                'descriptorsReference','descriptorsLocationReference',...
                'xQ','yQ','in','oceanMask','width','height','hReference');
        else
            load([resultFolder,'initial.mat'])
        end    
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

        % extract X & Y coordinates of the boundary of the changed area
        [xBoundaries,yBoundaries] = change_detection(currentImage,previousImage,segmentationMask);
        
        % measure the real size of each changed region
        [calvingSize,calvingSizePart] = measure_area2(xBoundaries,yBoundaries,xQ,yQ,hReference);
        
        % reformat the output to match the format of benchmark data
        calving_i{i} = calvingSize;
        coord_i{i} = cellfun(@horzcat,xBoundaries,yBoundaries,'UniformOutput',false);
        [~,imName,~] =  fileparts(imageFileName{i});
        calvName_i{i}=strjoin({imName,num2str(i),datestr(datenum(imageFileDate(i)),'ddmmyyyy_HHMM')},'_');
        calving_part_i{i} = calvingSizePart;
        disp(strcat(calvName_i{i},' done')); 
        lastit = i;
    end
end

%saving the data
save([resultFolder,'currentIm.mat'],'currentImage','currentMask','lastit');

calving_auto = cell2struct(calving_i,calvName_i,2);
calving_coord_auto = cell2struct(coord_i,calvName_i,2);
calvSizePart_auto = cell2struct(calving_part_i,calvName_i,2);
calving_cv_image_auto = cell2struct(cv_image,calvName_i,2);
calving_cv_both_auto = cell2struct(cv_both,calvName_i,2);
calving_currentImage_auto = cell2struct(currentImage_i,calvName_i,2);
calving_currentMask_auto = cell2struct(currentMask_i,calvName_i,2);
save([resultFolder,'calv_auto.mat'],'calving_auto');
save([resultFolder,'calvSizePart_auto.mat'],'calvSizePart_auto');
save([resultFolder,'calv_coord_auto.mat'],'calving_coord_auto');
save([resultFolder,'calv_current_auto.mat'],'calving_currentImage_auto','calving_currentMask_auto');
save([resultFolder,'calv_cv_auto.mat'],'calving_cv_image_auto','calving_cv_both_auto');

end