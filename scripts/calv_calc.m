clear all
close all
addpath('./Automatic/');
addpath('./Manual/');
addpath('./externalLibrary/');
addpath('./externalLibrary/efficientLBP/');
addpath('../config/');

% load config file. Config contains parameters for some methods
load config;
%define config as global variable
global config;

% Load calibration variables
%load calibration
save_calib = matfile(config.save_calib_mat, 'Writable', true);
    
% Init numbers    
load([config.resultdir,'/i_init.mat']);

% Screen Size
scrsz = get(groot,'ScreenSize');

% Save .mat file
save_calv = matfile(config.save_calv_mat, 'Writable', true);

if ~exist([config.resultdir,'currentIm.mat'])
    save_calv.pxSize_auto_i = zeros(length(config.imageFileName),config.part_num+1)*NaN;
    %save_calv.currentPoly_i = cell(length(config.imageFileName),1);
    save_calv.calving_part_i = zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calv.cv_image =  zeros(length(config.imageFileName),config.part_num+1);
    save_calv.mean_image =  zeros(length(config.imageFileName),config.part_num+1);
    save_calv.mean_both =  zeros(length(config.imageFileName),config.part_num+1);
    save_calv.cv_both = zeros(length(config.imageFileName),config.part_num+1);
    save_calv.ncalv_i = zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calv.calvLapse_i = zeros(length(config.imageFileName),1);
    

    calvingSize_i = cell(length(config.imageFileName),1);
    tForm_i = cell(length(config.imageFileName),1);
    coord_i = cell(length(config.imageFileName),1);
    
    it_init = i_init_mat{1};
    lastcalv = i_init_mat{1};
    lastit = i_init_mat{1}+1;
    it_sat = 1;
    nocalv = 0;

else
    load(config.save_coord);
    load(config.save_size);
    load(config.save_tForm);
    
    load([config.resultdir,'currentIm.mat']);
end

load([config.resultdir,'/initial',num2str(it_init),'.mat']);
% load satellite data reference
if config.calc_size == true
    frontfolder_it_sat = config.satdatefolder{it_sat};
else
    frontfolder_it_sat = '';
end
frontdir_it_sat = [config.frontdir,frontfolder_it_sat];
load ([frontdir_it_sat,'/H.mat'],'actualPixelSize');
actualPixelSize  = meshgrid(actualPixelSize,zeros(height,1));

for i=lastit:length(config.imageFileName)

    disp(strcat(config.calvName{i},' start'));
    % Check if there is an initiali.mat for next loop
    if ismember(i,cell2mat(i_init_mat))
      it_init = i;
    end
    
    im = str2double(regexp(config.imageFileName{i},'(\d*)','match'));

    previousImage = currentImage;
    previousMask = currentMask;    

    currentImage0 = imread(config.imagefiledir{i});
    currentImage0 = rgb2gray(currentImage0);
    currentImage = imcrop(currentImage0,config.cropRect);
        
    % Calculate the coefficient of variation
    [cv_prev_i,cv_cur_i,cv_both_i,mean_cur_i,mean_both_i] = cvIm(previousImage,currentImage,previousMask,width);

    if all(mean_cur_i(2:end) > save_calib.xmhigh(1,2:end)) || all(cv_cur_i(2:end) < save_calib.xcvlow(1,2:end))                    
        nocalv = i;                      
        tForm = tForm_i{i-1};
        save_calv.currentImage = previousImage;
        save_calv.currentMask = previousMask;
            
    else
        
        if config.nsat > 1 && it_sat < length(config.satdatefolder) && config.calv_size==true
            frontfolder_it_sat = config.satdatefolder{it_sat+1};
          frontdir_it_sat = [config.frontdir,frontfolder_it_sat];
	    Date_it = datetime(frontfolder_it_sat,'InputFormat','yyyy-MM-dd');
            if Date_it <= config.imageFileDate{i}
                it_sat = it_sat+1;
                % load satellite data reference
                load ([frontdir_it_sat,'/H.mat'],'actualPixelSize');
                actualPixelSize  = meshgrid(actualPixelSize,zeros(height,1));
            end
        end

        oceanMask0 = activecontour(currentImage0,initialoceanMask,config.activeContourIteration);
        oceanMask = imcrop(oceanMask0,config.cropRect);

        %register the current image using the first image as reference
        [currentImage0,tForm,prev] = registration4(descriptorsReference,descriptorsLocationReference,currentImage0,tForm_i{i-1});
        currentImage = imcrop(currentImage0,config.cropRect);
        if prev == 1
            currentMask = previousMask;
        else
            %extract SURF features, descriptors and descriptors' location of
            %the first image in a sequence
            surfReference = detectSURFFeatures(currentImage0);
            [descriptorsReference,descriptorsLocationReference] = extractFeatures(currentImage0,surfReference);

            % segment the current image based on the initial mask selected on
            % the first image.
            currentMask = activecontour(currentImage,initialMask,config.activeContourIteration);
            currentMask(oceanMask == 1) = 0;
            %currentMask(:,1:80) = 0;
        end
            [xcurrentPoly,ycurrentPoly] = mask2poly(currentMask,0,0);
            %save_calv.currentPoly_i{i} = horzcat(xcurrentPoly,ycurrentPoly);
       
        %combine the previous segmentation mask and current one to form
        %combined segmentation mask
        segmentationMask = previousMask & currentMask;

        % Calculate the coefficient of variation
        [cv_prev_i,cv_cur_i,cv_both_i,mean_cur_i,mean_both_i] = cvIm(previousImage,currentImage,segmentationMask,width);

        % extract X & Y coordinates of the boundary of the changed area
        if nocalv ~= i-1
            [xBoundaries,yBoundaries,calvMask_auto] = change_detection2(currentImage,previousImage,segmentationMask);
            pxSize_auto = zeros(1,config.part_num+1);
            for j=1:config.part_num
                pxSize_auto(j+1) = sum(sum(calvMask_auto(:,(j-1)*floor(width/config.part_num)+1:j*floor(width/config.part_num))));
            end
            calvMask_auto(segmentationMask == 0) = 99999;
            
            if all(pxSize_auto(2:end)>save_calib.maxcalv(1,2:end))
                nocalv = i;                      
                tForm = tForm_i{i-1};
                save_calv.currentImage = previousImage;
                save_calv.currentMask = previousMask;
            else
                        
                polys = cellfun(@horzcat,xBoundaries,yBoundaries,'UniformOutput',false);

                if length(polys(~cellfun('isempty',polys)))>0        
                    save_calv.calvLapse_i(i,1) = (datenum(config.imageFileDate{i})-datenum(config.imageFileDate{lastcalv}))*24*60;
                    lastcalv = i;

                    if isempty(polys)
                        coord_i{i} = {[]};
                        calvingSize_i{i} = {[]};
                        save_calv.calving_part_i(i,:) = zeros(1,config.part_num+1);
                        save_calv.ncalv_i(i,:) = zeros(1,config.part_num+1);
                        save_calv.pxSize_auto_i(i,:) = zeros(1,config.part_num+1);
                    else
                         [calvingSize_auto,calvingSizePart_auto,ncalv] = measure_area3(xBoundaries,yBoundaries,xQ,yQ,actualPixelSize);
                         calvingSize_i{i} = calvingSize_auto;
                         coord_i{i} = polys;
                         for p=2:config.part_num
                            if mean_cur_i(p) <= save_calib.xmhigh(1,p) && cv_cur_i(p) >= save_calib.xcvlow(1,p) && pxSize_auto(p) <= save_calib.maxcalv(1,p)
                                for k = 1:length(xBoundaries)
                                     in = poly2mask(xBoundaries{k},yBoundaries{k},height,width);%inpolygon(xQ,yQ,xBoundaries{k},yBoundaries{k});
                                     in_p = in(:,(p-1)*floor(width/config.part_num)+1:p*floor(width/config.part_num));
                                     if sum(in_p(:))<=config.alphaShapeRegionThreshold
                                        coord_i{j}{k} = [];
                                        calvingSize_i{j}{k} = [];
                                     else
                                        save_calv.pxSize_auto_i(i,p) = pxSize_auto(p);            
                                        % measure the real size of each changed region
                                        save_calv.calving_part_i(i,p) = calvingSizePart_auto(p);
                                        save_calv.ncalv_i(i,p) = ncalv(p);
                                     end
                                end
                            end
                         end
                    end
                else
                    coord_i{i} = {[]};
                    calvingSize_i{i} = {[]};
                    save_calv.calving_part_i(i,:) = zeros(1,config.part_num+1);
                    save_calv.ncalv_i(i,:) = zeros(1,config.part_num+1);
                    save_calv.pxSize_auto_i(i,:) = zeros(1,config.part_num+1);
                end
            end
        end

        lastit = i;
        tForm_i{i} = tForm;
        
        save(config.save_coord,'coord_i');
        save(config.save_size,'calvingSize_i');
        save(config.save_tForm,'tForm_i');

        save([config.resultdir,'currentIm.mat'],'currentImage','currentMask','lastit','it_init','actualPixelSize',...
        'it_sat','descriptorsReference','descriptorsLocationReference','lastcalv','nocalv');

       % Plot image and calving
        f = figure;
        set(f,'OuterPosition',[1 60 scrsz(3) scrsz(4)/2-60])        
        imshow(currentImage)
        hold on
        plot(xcurrentPoly,ycurrentPoly,'b')
        plot([xcurrentPoly(1),xcurrentPoly(end)],[ycurrentPoly(1),ycurrentPoly(end)],'b')
        hold on
        polys = coord_i{i};
        if ~isempty(polys)
            if length(polys{1})>0        
                for j =1:length(xBoundaries)
                    plot(xBoundaries{j},yBoundaries{j},'r','linewidth',2)
                end
            end
        end
        axis('equal')
         [~,imName,~] =  fileparts(config.imageFileName{i});
         set(gcf, 'PaperPosition', [0 0 30 10]); 
        print([config.figdir,datestr(datenum(config.imageFileDate{i}),'yyyymmdd_HHMM'),'_',imName],'-dpng','-r0')
        close(f)
    end
    
    % Save image characteristics
    save_calv.cv_image(i,:) =  cv_cur_i;
    save_calv.mean_image(i,:) =  mean_cur_i;
    save_calv.cv_both(i,:) = cv_both_i;
    save_calv.mean_both(i,:) =  mean_both_i;

    % Check if there is an initiali.mat for next loop
    if ismember(i,cell2mat(i_init_mat))
       load([config.resultdir,'/initial',num2str(i),'.mat']);
    end
    disp(strcat(config.calvName{i},' done'));
end
