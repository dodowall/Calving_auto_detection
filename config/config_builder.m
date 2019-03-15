config = struct;

config.maindir = fileparts(fileparts(mfilename('fullpath')));

%%
%%%%%%%%%%%%%%%%%%%%%
% To be changed     %
%%%%%%%%%%%%%%%%%%%%%
% Name of project
config.projname = '2014';

% Folder with satellite images to calculate distance from the front
config.satdir = [config.maindir,'/',config.projname,'/landsat/'];

% true if you want to calculate the size, false otherwise
config.calc_size = false;

%% Image and camera parameters
%Image pixel size
config.imagePixelSize = 5e-6;
%Focal length of the lens
config.focalLength = 28e-3;
% Camera sensor size
config.xsens = 22.2;
config.ysens = 14.8;
% Crop image around the front
config.cropRect = [169 1189 3911 443];
% Number of zones
config.part_num = 6;
% Time lapse (min)
config.tl = 14;

%% Camera and satellite image parameters
% Zoom in the satellite image corners
config.ximz1 = 551530;
config.ximz2 = 554450;
config.yimz1 = 8707800;
config.yimz2 = 8711000;

% Position of front on the image
config.xl1 = 552080;
config.xl2 = 554000;
config.yl1 = 8710140;
config.yl2 = 8708100;
 
% Camera position
config.xcam = 550686;
config.ycam = 8709729;
config.zcam = 73;

% Limits of camera view
config.xlimref1 = 560910;
config.ylimref1 = 8704700; 
%% Parameters for the automatic detection
% the default size of maximum calving (the image is not taken into account
% above this value)
config.max_calv_default = 3000;
% Define how often the contour must be manually set (number of pics in
% between)
config.contour_init = 1000; % depends on the time lapse and the glacier velocity

%Filter size and radius for Local Binary Pattern
config.LBPFilterSize = 20;%8;
config.LBPFilterRadius = 5;%5;

%Filter size for mean and median filter used for merging two type of
%difference image
config.meanFilterSize = [11 11];
config.medianFilterSize = [3 3];

%Weight for merging two difference image
config.mergeAlpha = 0.1;

%Window size and constant for Adaptive threshold
config.adaptWindow = 25;
config.adaptConstant = 0.08;
%Type of adaptive threshold, 0 for mean, 1 for median
config.adaptFilterType = 1;
%Radius and minimum region area of Alpha Shape
config.alphaShapeRadius = 10;
config.alphaShapeRegionThreshold = 100;

%Number of iterations for Active Contour
config.activeContourIteration = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Folders and files                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each satellite image should be in a folder with the date
satdatefolder = dir(config.satdir);
satdatefolder = {satdatefolder.name}';
config.satdatefolder = satdatefolder(3:end);
config.nsat = length(config.satdatefolder);

% Folder with camera images
config.imagedir = [config.maindir,'/',config.projname,'/images/'];
imfold = dir(fullfile(config.imagedir));
config.imagefolders = sort({imfold(3:end).name}');
%image format and prefix
config.imagePrefix = 'IMG';
config.imageFormat = 'JPG';

% List of characteristics
config.fogList = {'Normal','Dense fog','High illumination'};

% Results folder
config.resultdir = [config.maindir,'/results_',config.projname,'/'];

% Saved variables
config.calibdir = [config.resultdir,'calibration/'];
if ~exist(config.calibdir)
    mkdir(config.calibdir);
end
config.save_calib_mat = [config.calibdir,'save_calib_',config.projname,'.mat'];
config.save_calib_tForm = [config.calibdir,'save_calib_tForm_',config.projname,'.mat'];

config.save_calv_mat = [config.resultdir,'save_calv_',config.projname,'.mat'];
config.save_coord = [config.resultdir,'save_coord_',config.projname,'.mat'];
config.save_tForm = [config.resultdir,'save_tForm_',config.projname,'.mat'];
config.save_size = [config.resultdir,'save_size_',config.projname,'.mat'];



% Folder where front data are stored after running front_position
config.frontdir = [config.resultdir,'/front/'];
if ~exist(config.frontdir)
    mkdir(config.frontdir);
end

for i=1:length(config.satdatefolder)
    config.frontdir_i = [config.frontdir,config.satdatefolder{i}];
    if ~exist(config.frontdir_i)
        mkdir(config.frontdir_i);
    end
end

% Folder with resulting figures
config.figdir = [config.resultdir,'/figures/'];
if ~exist(config.figdir)
    mkdir(config.figdir);
end

config.imagefiledir = [];
config.imageFileDate = [];
config.imageFileName = [];
config.calvName = [];
for im=1:length(config.imagefolders)
    imPath = [config.imagedir,config.imagefolders{im},'/'];
    imageFiles = dir(fullfile(imPath,strcat('*.',config.imageFormat)));
    imFileName = sort({imageFiles.name}');
    imageName = cell(length(imFileName),1);
    for i=1:length(imFileName)
        [~,imName,~] =  fileparts(imFileName{i});
        imageName{i} = imName;
    end
    config.imagefiledir = [config.imagefiledir;cellfun(@(x) [imPath,'/',x,'.',config.imageFormat],imageName,'UniformOutput',false)];
    config.imageFileName = [config.imageFileName;cellfun(@(x) [config.imagefolders{im},'_',x],imageName,'UniformOutput',false)];

    for i=1:length(imFileName)
        tmp = imfinfo([imPath,imFileName{i}]);
        imFileDate = datetime(datevec(tmp.DateTime,'yyyy:mm:dd HH:MM:SS'));
        config.imageFileDate = [config.imageFileDate;cellfun(@(x) x,{imFileDate},'UniformOutput',false)];  
        calvName_i=strjoin({imageName{i},num2str(i),datestr(datenum(imFileDate),'ddmmyyyy_HHMM')},'_');
        config.calvName = [config.calvName;cellfun(@(x) x,{calvName_i},'UniformOutput',false)];        
    end
end

save([config.maindir,'/config/config'],'config');
