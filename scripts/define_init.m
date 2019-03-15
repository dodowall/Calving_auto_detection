clear all
close all
addpath('./Automatic/');
addpath('./Manual/');
addpath('./externalLibrary/');
addpath('../config/');

% load config file. Config contains parameters for some methods
load config;
%define config as global variable
global config;

% Check if there is an initiali.mat
initfiles = dir(fullfile(config.resultdir,'initial*.mat'));
if isempty(initfiles)
    it = 1;
    i_init_mat = cell(floor(length(config.imageFileName)/config.contour_init)+1,1);
    if config.calc_size == true
        % Load front/ocean contour
        if config.calc_size == true
            frontfolder_it_sat = config.satdatefolder{1};
        else
            frontfolder_it_sat = '';
        end
        load ([frontfolder_it_sat,'/front.mat'],'ximp','yimp');
        x_ocean = ximp+config.cropRect(1);
        y_ocean = yimp+config.cropRect(2);
        pts_ocean = 45;
    end
else
    initfilesname = {initfiles.name}';
    initnum = zeros(length(initfilesname),1);
    for i=1:length(initfilesname)
        [~,initname,~] =  fileparts(initfilesname{i});
        initnum(i) = str2num(erase(initname,'initial'));
    end
    it = max(initnum);
    
    load([config.resultdir,'/initial',num2str(it),'.mat'],'initialoceanMask','initialMask');
    [x_ocean,y_ocean] = mask2poly(initialoceanMask,0,0);
    [x,y] = mask2poly(initialMask,0,0);
    pts_ocean = 30;
    it =  max(initnum) + config.contour_init;
    load([config.resultdir,'/i_init.mat']);
end

% Screen Size
scrsz = get(groot,'ScreenSize');

j = 1;
for i=it:config.contour_init:length(config.imageFileName)
    info = helpdlg(...
         {'Choose the image that is best to define initial contour.',...
         'Then draw the ocean contour.',...
         'Check and press enter.',...
         'Then draw the front contour.'...
         'Check and press enter.',...
         'Next image!'},'Info');

     uiwait(info);
 
    % Choose image that is good to start with
    i_init = i;
    closef = false;
    while closef == false
        f = figure('OuterPosition',[1 60 scrsz(3) scrsz(4)-60],'Toolbar','figure','Name',config.imageFileName{i_init});
        h1 = uicontrol('Style', 'pushbutton', 'String', 'Previous image',...
        'FontSize',20,'Position', [100 2*scrsz(4)/3 200 100],...
        'Callback', 'i_init=max(1,i_init-1);uiresume(gcbf);');            
        h2 = uicontrol('Style', 'pushbutton', 'String', 'Next image',...
        'FontSize',20,'Position', [100 scrsz(4)/3 200 100],...
        'Callback', 'i_init=i_init+1;uiresume(gcbf);');        
        h3 = uicontrol('Style', 'pushbutton', 'String', 'Take this one!',...
            'Position', [100 100 200 100],...
            'FontSize',20,'Callback', 'closef=true;uiresume(gcbf);'); 

        im = str2double(regexp(config.imageFileName{i_init},'(\d*)','match'));
        %read image and convert to grayscale
        currentImage0 = imread(config.imagefiledir{i_init});
        currentImage0 = rgb2gray(currentImage0);

        imshow(currentImage0);

        uiwait(gcf);
        close(f);
        
    end

    [height0,width0]=size(currentImage0);

    %extract SURF features, descriptors and descriptors' location of
    %the first image in a sequence
    surfReference = detectSURFFeatures(currentImage0);
    [descriptorsReference,descriptorsLocationReference] = extractFeatures(currentImage0,surfReference);

    %show the image and let user select ocean
    figure('Name','Draw the ocean contour and enter');
    imshow(currentImage0);
    hold on
    %[xocean,yocean]=getline;
    if config.calc_size == true
        [xim_1,yim_1] = interp_pts(x_ocean,y_ocean,pts_ocean);
        linetoedit = line(xim_1,yim_1);
        [x_ocean,y_ocean] = editpolyline(linetoedit);
    else
        [x_ocean,y_ocean]=getline;
    end
    close;
    
    initialoceanMask = poly2mask(x_ocean,y_ocean,height0,width0);
    oceanMask0 = activecontour(currentImage0,initialoceanMask,config.activeContourIteration);
    [x_ocean,y_ocean] = mask2poly(oceanMask0,0,0);        

    f = figure('Name','Check and close');
    imshow(currentImage0)
    hold on
    plot(x_ocean,y_ocean)
    uiwait(gcf);
    close;

    %Crop image, the cropped size is matched the benchmark image size
    currentImage = imcrop(currentImage0,config.cropRect);
    oceanMask = imcrop(oceanMask0,config.cropRect);
    [height,width]=size(currentImage);

    [xQ,yQ]=meshgrid(1:width,1:height);

    %show the image and let user select an initial segmentation mask
    figure('Name','Draw the front contour and enter');
    imshow(currentImage);
    if it == 1
        [x,y]=getline;
    else
        [x1,y1] = interp_pts(x,y,30);
        linetoedit = line(x1,y1);
        [x,y] = editpolyline(linetoedit);
    end
    close;
        
    %segment the image usign Chan-vese algorithm based on the initial
    %mask
    initialMask = poly2mask(x,y,height,width);
    currentMask = activecontour(currentImage,initialMask,config.activeContourIteration);
    % Remove points situated in the ocean
    currentMask(oceanMask == 1) = 0;

    % Points situated at the front and in the ocean
    [x_poly,y_poly] = mask2poly(currentMask,config.cropRect(1), config.cropRect(2));        

    figure('Name','Check and close');
    imshow(currentImage0)
    hold on
    plot(x_poly,y_poly)
    plot(x_ocean,y_ocean)
    uiwait(gcf);
    close;

    currentMask0 = poly2mask(x_poly,y_poly,height0,width0);
    mask0 = currentMask0 | oceanMask0;
    [x0,y0] = mask2poly(mask0,0,0);          
    in = inpolygon(descriptorsLocationReference.Location(:,1),descriptorsLocationReference.Location(:,2),x0,y0);

    save([config.resultdir,'/initial',num2str(i_init),'.mat'],'initialMask','currentImage',...
        'currentMask','initialoceanMask',...
        'descriptorsReference','descriptorsLocationReference',...
        'xQ','yQ','in','oceanMask','width','height');
    
   i_init_mat{j} = i_init;
   j = j+1;
end  
save([config.resultdir,'/i_init.mat'],'i_init_mat');
