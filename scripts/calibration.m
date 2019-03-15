clear all
close all
addpath('./Automatic/');
addpath('./Manual/');
addpath('./externalLibrary/');
addpath('./externalLibrary/efficientLBP/');
addpath('../config/');

choice_list = {'YES','NO'};

% load config file. Config contains parameters for some methods
load config;
%define config as global variable
global config;

% Init numbers    
load([config.resultdir,'/i_init.mat']);

% Screen Size
scrsz = get(groot,'ScreenSize');

% Save .mat file
save_calib = matfile(config.save_calib_mat, 'Writable', true);

if ~exist([config.calibdir,'currentCalib.mat'])
    save_calib.cv_image =  zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calib.mean_image =  zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calib.ncalv_i = zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calib.calving_part_i =  zeros(length(config.imageFileName),config.part_num+1)*NaN;
    save_calib.pxSize_auto_i = zeros(length(config.imageFileName),config.part_num+1)*NaN;

    tForm_i = cell(length(config.imageFileName),1);
    
    save_calib.fogPart_i =  ones(1,config.part_num+1);
    save_calib.choice_i = ones(1,config.part_num+1);
    
    save_calib.xmhigh = zeros(1,config.part_num+1)*NaN;
    save_calib.xcvlow = zeros(1,config.part_num+1)*NaN;
    save_calib.xmhigh_std = zeros(1,config.part_num+1)*NaN;
    save_calib.xcvlow_std = zeros(1,config.part_num+1)*NaN;
    save_calib.maxcalv = ones(1,config.part_num+1)*config.max_calv_default;
    
    it_init = i_init_mat{1};
    lastcalv = i_init_mat{1};
    lastit = i_init_mat{1}+1;
    it_sat = 1;
    imNumtostart = i_init_mat{1}+1;

else
    load(config.save_calib_tForm);
    load([config.calibdir,'currentCalib.mat']);
    
    % Construct a questdlg with three options
    choice = questdlg('Where do you want to start again?', 'Start menu',...
	'Where I ended','Choose...','Where I ended');

    % Handle response
    switch choice
        case 'Where I ended'
            imNumtostart = lastit;
            
        case 'Choose...'
            choice2 = questdlg('What to choose?', 'Choice',...
	'From image number','From image indice','From image date','From image number');
        switch choice2
            case 'From image number'
                which_num = inputdlg('Image number','Image number to start with',[1 50],config.imageFileName(min(lastit,length(config.imageFileName))));        
                imNumtostart = find(strcmp(config.imageFileName,strjoin(['IMG_',which_num,'.JPG'],'')));

            case 'From image indice'
                which_ind = inputdlg('Image number','Image indice to start with',[1 50],{num2str(lastit)});        
                imNumtostart = str2double(which_ind{1});

            case 'From image date'

                DateVector = datevec(imageDate(min(lastit,length(config.imageFileName))),'dd-mmm-yyyy HH:MM:ss');             
                defaultanswer = {num2str(config.imageFileDate{1}),num2str(config.imageFileDate{2}),num2str(config.imageFileDate{3}),num2str(config.imageFileDate{4}),num2str(config.imageFileDate{5}),num2str(config.imageFileDate{6})};
                which_date = inputdlg({'Year (yyyy)','Month (mm)','Day (dd)','Hour (hh)','Minute (MM)','Seconds (ss)'},'Image date to start with',[1 50;1 50;1 50;1 50;1 50;1 50],defaultanswer);
                datetostart = datestr([str2double(which_date{1}),str2double(which_date{2}),str2double(which_date{3}),str2double(which_date{4}),str2double(which_date{5}),str2double(which_date{6})]); 
                imNumtostart = find(strcmpi(config.imageFileDate,datetostart));
        end
    end
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

nocalv = 0;
it = 1;
it_sat = 1;

% Part names
for p=1:config.part_num
    partNames{p} = ['Part ',num2str(p)];
end

%%
%%%%%%%%%%%%%%%%%%
% Start picking  %
%%%%%%%%%%%%%%%%%%
i=imNumtostart;
global im
global fogPart;
global choice_im;
im = 0;

% Number of images
n = length(config.imageFileName);

while i<=n
    disp(strcat(config.calvName{i},' start'));
    % Check if there is an initiali.mat for next loop
    if ismember(i,cell2mat(i_init_mat))
      it_init = i;
    end
    
    im1 = str2double(regexp(config.imageFileName{i},'(\d*)','match'));

    previousImage = currentImage;
    previousMask = currentMask;    

    currentImage0 = imread(config.imagefiledir{i});
    currentImage0 = rgb2gray(currentImage0);
    currentImage = imcrop(currentImage0,config.cropRect);

    oceanMask0 = activecontour(currentImage0,initialoceanMask,config.activeContourIteration);
    oceanMask = imcrop(oceanMask0,config.cropRect);

%     %register the current image using the first image as reference
%     [currentImage0,tForm,prev] = registration4(descriptorsReference,descriptorsLocationReference,currentImage0,tForm_i{i-1});
%     currentImage = imcrop(currentImage0,config.cropRect);
%     if prev == 1
%         currentMask = previousMask;
%     else
%         %extract SURF features, descriptors and descriptors' location of
%         %the first image in a sequence
%         surfReference = detectSURFFeatures(currentImage0);
%         [descriptorsReference,descriptorsLocationReference] = extractFeatures(currentImage0,surfReference);

        % segment the current image based on the initial mask selected on
        % the first image.
        currentMask = activecontour(currentImage,initialMask,config.activeContourIteration);
        currentMask(oceanMask == 1) = 0;
        %currentMask(:,1:80) = 0;
%    end
        [xcurrentPoly,ycurrentPoly] = mask2poly(currentMask,0,0);
        %save_calib.currentPoly_i{i} = horzcat(xcurrentPoly,ycurrentPoly);

    %combine the previous segmentation mask and current one to form
    %combined segmentation mask
    segmentationMask = previousMask & currentMask;

    % Calculate the coefficient of variation
    [cv_prev_i,cv_cur_i,cv_both_i,mean_cur_i,mean_both_i] = cvIm(previousImage,currentImage,segmentationMask,width);
    % Save image characteristics
    save_calib.cv_image(i,:) =  cv_cur_i;
    save_calib.mean_image(i,:) =  mean_cur_i;
    
    a=size(currentImage);
  %  part_size = a(2)/config.part_num;
  %  x_lim = [0.5 a(2)];
    title_im = config.calvName{i};
        
    %%
    
    % Decide exposition / fog 
%     if i==imNumtostart
%         info = helpdlg(...
%              {''},'Info');
%          uiwait(info);
%     end
    if sum(~isnan(save_calib.fogPart_i(:,1)))<i
        fogPart = save_calib.fogPart_i(i-1,:);
        choice_im = save_calib.choice_i(i-1,:);
    else
        fogPart = save_calib.fogPart_i(i,:);
        choice_im = save_calib.choice_i(i,:);
        im = 1;
    end
    f=figure('OuterPosition',[1 60 scrsz(3) scrsz(4)-60],'Toolbar','figure','Name',title_im);
    handles = guidata(gcf);
    for p=1:config.part_num         
        if i == i_init_mat{1}+1
            fog_i = ones(1,config.part_num+1);
            ch_i = ones(1,config.part_num+1);
        else
            fog_i = save_calib.fogPart_i(i-1+im,p);
            ch_i = save_calib.choice_i(i-1+im,p);
        end

        handles.bg{p} = uibuttongroup('Title',partNames{p},...
                  'Position',[1/config.part_num*(p-1) 0.85 config.part_num 0.15],...
                  'SelectionChangedFcn',{@getFog,p,config.fogList,config.part_num});                               
        
        for fl=1:length(config.fogList)
             handles.hfl{p,fl} = uicontrol(handles.bg{p},'Style', 'radiobutton', 'String', config.fogList{fl},...
            'Position', [(scrsz(3))/config.part_num^2 fl*20 100 20]); 
            guidata(gcf,handles)  
        end
        handles.bg{p}.SelectedObject = handles.hfl{p,fog_i}; 
        
        handles.bg2{p} = uibuttongroup('Title','Choice',...
                  'Position',[1/config.part_num*(p-1) 0.75 config.part_num 0.1],...
                  'SelectionChangedFcn',{@getChoice,p,choice_list,config.part_num});                               
        
        for ch=1:length(choice_list)
             handles.hfl2{p,ch} = uicontrol(handles.bg2{p},'Style', 'radiobutton', 'String', choice_list{ch},...
            'Position', [(scrsz(3))/config.part_num^2 ch*20 100 20]); 
            guidata(gcf,handles)  
        end
        handles.bg2{p}.SelectedObject = handles.hfl2{p,ch_i}; 
        
        guidata(gcf,handles)  
    end   
    h1 = uicontrol('Style', 'pushbutton', 'String', 'Previous image',...
        'Position', [1*scrsz(3)/6 80 100 20],...
        'Callback', 'im=1;uiresume(gcbf);');   
    h2 = uicontrol('Style', 'pushbutton', 'String', 'Next image',...
        'Position', [1*scrsz(3)/2 80 100 20],...
        'Callback', 'im=0;uiresume(gcbf);');       
    h3 = uicontrol('Style', 'pushbutton', 'String', 'Stop iterations',...
        'Position', [5*scrsz(3)/6 80 100 20],...
        'Callback', 'im=-1;uiresume(gcbf);'); 
    set(f,'WindowKeyPressFcn',{@fog_pressed_fcn,config.fogList,config.part_num});
    
    % First plot
    subplot('Position',[0.13 0.42 0.775 0.3])
    imshow(currentImage);
    for k=1:config.part_num
        rectangle('Position',[(k-1)*(a(2))/config.part_num,0.5,(a(2))/config.part_num,a(1)]);
%         text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5,['Kappa = ',num2str(kappa_i(i,k+1))]);
%         text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5+75,['MCC = ',num2str(mcc_i(i,k+1))]);
%         text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5+150,['AreaDiff = ',...
%             num2str((pxSize_auto_i(i,k+1)-pxSize_manual_i(i,k+1)))]);
%         text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5+225,['SizeMeasure = ',...
%             num2str(size_i(i,k+1))]);
        text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5,['cv = ',num2str(cv_cur_i(k+1))]);
        text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5+75,['mean = ',num2str(mean_cur_i(k+1))]);
        text((k-1)*(a(2))/config.part_num,scrsz(4)/1.5+150,['std = ',num2str(cv_cur_i(k+1)*mean_cur_i(k+1))]);
   end
    hold on
    % Extract poly manual
%     polys = coord_manual_i{i};
% 
%     if length(polys{1})>0
%         xBoundaries = cellfun(@(x) x(:,1),polys, 'UniformOutput', false);
%         yBoundaries = cellfun(@(x) x(:,2),polys, 'UniformOutput', false);
%         xBoundaries = xBoundaries(~cellfun('isempty',xBoundaries));       
%         yBoundaries = yBoundaries(~cellfun('isempty',yBoundaries));
%  
%         for j =1:length(xBoundaries)
%             plot(xBoundaries{j},yBoundaries{j},'linewidth',2)
%         end
%     end
   
    hold on 
    
    test_xmhigh = save_calib.xmhigh-save_calib.xmhigh_std;
    test_xcvlow = save_calib.xcvlow+save_calib.xcvlow_std;
    
    if any(mean_cur_i>=test_xmhigh) | any(cv_cur_i<=test_xcvlow) | i==i_init_mat{1}+1

        subplot('Position',[0.13 0.11 0.775 0.3])
        imshow(currentImage);
        for k=1:config.part_num
            rectangle('Position',[(k-1)*(a(2))/config.part_num,0.5,(a(2))/config.part_num,a(1)]);
        end
        hold on
                
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
        %save_calib.currentPoly_i{i} = horzcat(xcurrentPoly,ycurrentPoly);

        %combine the previous segmentation mask and current one to form
        %combined segmentation mask
        segmentationMask = previousMask & currentMask;

        if config.nsat > 1 && it_sat < length(config.satdatefolder) && config.calc_size==true
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

        % extract X & Y coordinates of the boundary of the changed area
        if nocalv ~= i-1
            [xBoundaries,yBoundaries,calvMask_auto] = change_detection2(currentImage,previousImage,segmentationMask);
            pxSize_auto = zeros(1,config.part_num+1);
            for j=1:config.part_num
                pxSize_auto(j+1) = sum(sum(calvMask_auto(:,(j-1)*floor(width/config.part_num)+1:j*floor(width/config.part_num))));
            end
            calvMask_auto(segmentationMask == 0) = 99999;

            polys = cellfun(@horzcat,xBoundaries,yBoundaries,'UniformOutput',false);

            if length(polys(~cellfun('isempty',polys)))>0
                if isempty(polys)
                    save_calib.calving_part_i(i,:) = zeros(1,config.part_num+1);
                    save_calib.ncalv_i(i,:) = zeros(1,config.part_num+1);
                    save_calib.pxSize_auto_i(i,:) = zeros(1,config.part_num+1);
                else
                    save_calib.pxSize_auto_i(i,:) = pxSize_auto;            
                    % measure the real size of each changed region
                    [calvingSize_auto,calvingSizePart_auto,ncalv] = measure_area3(xBoundaries,yBoundaries,xQ,yQ,actualPixelSize);
                    save_calib.calving_part_i(i,:) = calvingSizePart_auto';
                    save_calib.ncalv_i(i,:) = ncalv';  

                    for j =1:length(xBoundaries)
                        plot(xBoundaries{j},yBoundaries{j},'linewidth',2)
                    end

                end
            else
                coord_i{i} = {[]};
                calvingSize_i{i} = {[]};
                save_calib.calving_part_i(i,:) = zeros(1,config.part_num+1);
                save_calib.ncalv_i(i,:) = zeros(1,config.part_num+1);
                save_calib.pxSize_auto_i(i,:) = zeros(1,config.part_num+1);
            end
        end
    end

    
    %%
    uiwait(gcf);
    close(f);
    save_calib.fogPart_i(i,:) = fogPart;
    save_calib.choice_i(i,:) = choice_im;
    
    if any(save_calib.choice_i(i,:) == 2)
        tForm = tForm_i{i-1};
        currentImage = previousImage;
        currentMask = previousMask;
        nocalv = i;

    else
        lastcalv = i;
    end
    lastit = i;
    tForm_i{i} = tForm;

    save(config.save_calib_tForm,'tForm_i');
    
    if i > i_init_mat{1}+1 & any(save_calib.choice_i(:,p)==2)
        mean_image = save_calib.mean_image;
        cv_image = save_calib.cv_image;
        pxSize_auto = save_calib.pxSize_auto_i;
        for p=1:max(config.part_num)
            save_calib.xmhigh(:,p) = min(mean_image(save_calib.choice_i(:,p)==2 & save_calib.fogPart_i(:,p)==3,p));
            save_calib.xcvlow(:,p) = max(cv_image(save_calib.choice_i(:,p)==2 & save_calib.fogPart_i(:,p)==2,p));
            save_calib.xmhigh_std(:,p) = std(mean_image(save_calib.choice_i(:,p)==2 & save_calib.fogPart_i(:,p)==3,p));
            save_calib.xcvlow_std(:,p) = std(cv_image(save_calib.choice_i(:,p)==2 & save_calib.fogPart_i(:,p)==2,p));
            save_calib.maxcalv(:,p) = min(pxSize_auto(save_calib.choice_i(:,p)==2 & save_calib.fogPart_i(:,p)==1,p));
        end
    end
    
    save([config.calibdir,'currentCalib.mat'],'currentImage','currentMask','lastit','it_init','actualPixelSize',...
    'it_sat','descriptorsReference','descriptorsLocationReference','lastcalv','nocalv','tForm');

    % Check if there is an initiali.mat for next loop
    if ismember(i,cell2mat(i_init_mat))
       load([config.resultdir,'/initial',num2str(i),'.mat']);
    end
    disp(strcat(config.calvName{i},' done'));
        
    if im == 0
        i=i+1;
    elseif im == 1
        i = i-1;
    else
        i = n+1;
    end

end
