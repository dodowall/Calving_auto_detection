# Calving_detection
Program to automatically detect calving at the front of a tidewater glacier from timelapse images

how to download: https://zenodo.org/badge/latestdoi/175037496

how to cite: Vallot D.,  Adinugroho S., Strand R. and Pettersson R.: Program to automatically detect calving at the front of a tidewater glacier from timelapse images, Zenodo, doi: 10.5281/zenodo.2591197, 2019.

More details in the scientific article:


REQUIREMENT : 
The application uses some internal functions introduced in MATLAB 2014b. Therefore, it won’t work on MATLAB version lower than 2014b. 
Some functions are also provided by external libraries. For simplicity, I have included all external libraries in a folder named “externalLibrary”. Make sure you include the “externalLibrary” folder into your MATLAB path. It also requires the Deep Learning Toolbox.


*******************
* Input           
*******************
The main path should contain a folder with input data (images and satellite data if size calculation is wanted) called Tuna14 for example (the user is deciding the name that will be written in the config builder).
The images should be placed in Tuna14/images in one or several folders (for example input/images/CANON100, Tuna14/images/CANON200) but not directly in Tuna14/images.
The satellite data should be placed in Tuna14/landsat in one or several folders with date (for example Tuna14/landsat/2014-06-20). The name of the folder landsat can be changed in the config builder. In the example, it is just an xml file but the landsat data can be downloaded from EarthExplorer (https://earthexplorer.usgs.gov) for example.

*******************
* Run the scripts 
*******************
There are different steps to go through in order to automatically detect calving and (if wanted) calculate the size of calving. The four first steps are preparation steps and the last one do the job.
- 1st step: configuration
  open config/config_builder.m
  change parameters given your file structure, camera, images and satellite characteristics as well as the parameters for the automatic detection.
  ignore the satellite characteristics if no size calculation.
  run config/config_builder.m

- 2nd step: front position
  run scripts/front_postion.m even if you don't need the size calculation
  if you calculate the size, you will need to delineate the limit between the front and the ocean on the image and on the satellite image.

- 3rd step: initialisation
  run scripts/define_init.m
  depending on config.contour_init (defined in 1st step) and the number of images, you will have to manually delineate several initial fronts.

- 4th step: calibration
  run scripts/calibration.m
  you will go through each image and for each zone decide if it is normal, foggy or too illuminated and if you want to keep it or not
  the lower image shows the automatic detection result and if it is not realistic, don't keep the zone (it won't always appear depending on the characteristics of the image)

- 5th step: automatic detection
  run calv_calc.m

*******************
* Results         
*******************
For each image, the calving detection is plotted and stored in results_Tuna14/figures.
The coordinates of each calving event for each image are stored in results_Tuna14/save_coord_Tuna14.mat.
The size of each calving event for each image is stored in results_Tuna14/save_size_Tuna14.mat.
Other interesting characteristics of the calving detection per zone are stored in results_Tuna14/save_calv_Tuna14.mat.

