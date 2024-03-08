# Changelog


## 2024-03-08

 Version 2.3  
 **NB: this version changes the analysis methods and the results are different from previous versions**  
- Add fascicle detection outline at the end of analysis  
- Add possibility to let user test and adjust parameters for aponeurosis detection  
- Add manual cropping for panoramic scans  
- Revert fix for aponeuroses display in version 2.2.1 because it prevented analysis of panoramic scans  
- Improve aponeuroses and fascicle prefiltering by using MorphoLibJ plugin instead of custom directional filter  
- Improve fascicle prefiltering  
The sigma value of the vesselness filter now matches the chosen sigma used for the Laplacian of Gaussian. This make the prefiltering somewhat adapted to expected fascicle thickness.  
The "test" function for the Laplacian of Gaussian value was removed as it could not be maintained with this change.  
One of the filtering steps was also removed as it was not improving detection.  
- Fix various bugs  
	
## 2023-05-14

 Version 2.2.1  
- fix manual cropping  
- fix aponeuroses display in cases of relative misalignment  
- change settings interface  
	secondary settings are now shown in a secondary interface if necessary  
	settings are stored in a 'SMA_prefs.txt' file in Fiji 'plugins' directory  

## 2023-04-07

 Version 2.2  
  **NB: this version changes the analysis methods and the results are different from previous versions**  
- implement contrast enhancement (CLAHE plugin) for fascicles  
- implement handling for movies and image sequences  
- fix bug causing fascicle insertions to not be displayed correctly when the muscle is not horizontal  
 - various bug fixes  

## 2023-03-05

 Version 2.1  
  **NB: this version changes the analysis methods and the results are different from previous versions**  
 - implement option to open single image  
 - improve handling of DICOM metadata  
 - change vesselness function to detect fascicles after FFT from Tubeness to Frangi  
 - make aponeurosis enhancing step (CLAHE plugin) optional.   
 	This filter only helps in case of faint aponeuroses but could make the analysis fail in other cases  
 - add expected aponeurosis length option.  
 	relative to FoV length (default 80%), can be helpful in cases of false positives  
 - complete commenting  
 - add "developper" mode   
 	Disable batch mode and show intermediate steps of analysis in single file mode  
 
 Version 2.0.0  
 **NB: this version changes the analysis methods and the results are different from previous versions**  
 - add possibility to take fascicle curvature into account  
 	this is implemented by detecting fascicle orientation in superficial and deep regions and 
 	reconstructing a composite fascicle by simple spline fitting of two fascicle segments or 
 	by circle fitting. These methods are based on fitting curves across the insertion points and the intersection points between fascicles segments at 	   mid-distance between aponeuroses.  
 - detect fascicle orientation based on single (straight fascicle method) or paired  
 	(curved fascicles methods) ROIs. This change is possible because the image is now rotated 
 	by an angle equal to lower aponeurosis angle and ROIs can be broader.  
 - add panoramic mode  
 	currently, add the possibility to standardise ROI width by selecting a subregion with a 
 	width = FoV width * 'x field of view' factor, chosen by user.  
 	**NB: there is no detection of the original FoV width, a width of 5 cm is currently assumed!**  
 - various changes to filtering method before measurement of fascicle orientation.  
 - various bug fixes and code improvements.  

## 2022-05-04

Version 1.7.1  
- add possibility to extrapolate aponeuroses from 100% or 50% of detected length.  
- add fascicle length calculated from mean thickness and angle ("Fascicle_length_trig").   
	NB: The aponeurosis angle is the slope of the straight line fitted to the chosen portion of the aponeurosis (i.e. 50 or 100%). **It will influence 	   the obtained pennation angle and fascicle length with any method**.

## 2020-03-17

Version 1.7  
- add error handling.  
- check depencies.  
- fix bug that caused analysis to fail when images are flipped and cropping manual.  
