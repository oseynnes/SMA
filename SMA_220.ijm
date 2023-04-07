///////////////////////////////////////////////
// SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS //
///////////////////////////////////////////////
/*
*	Requires OrientationJ plugin v2.0.3 (http://bigwww.epfl.ch/demo/orientation - Z. P sp ki, M. Storath, D. Sage, M. Unser, "Transforms and Operators for Directional Bioimage Analysis: A Survey," Advances in Anatomy, Embryology and Cell Biology, vol. 219, Focus on Bio-Image Informatics, Springer International Publishing, May 21, 2016.)
*	Requires Canny Edge Detector plugin (https://imagej.nih.gov/ij/plugins/canny/index.html - Tom Gibara)
*	Requires Non Local Means Denoise plugin v1.4.6 (https://imagej.net/Non_Local_Means_Denoise - Pascal Behnel, Thorsten Wagner)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

////////////////
// Change log //
////////////////
/*
 * Version 2.2
 *  - implement contrast enhancement (CLAHE plugin) for fascicles
 *  - implement handling for movies and image sequences
 *  - fix bug causing fascicle insertions to not be displayed correctly when the muscle is not horizontal
 *  - various bug fixes
 * 
 * 	Version 2.1
 * 	- implement option to open single imnage
 * 	- improve handling of DICOM metadata
 * 	- change vesselness function to detect fascicles after FFT from Tubeness to Frangi
 * 		This is important as it means the results are different from previous versions
 * 	- make aponeurosis enhancing step (CLAHE plugin) optional.
 * 		This filter only helps in case of faint aponeuroses but could make the analysis fail in other cases
 * 	- add expected aponeurosis length option.
 * 		relative to FoV length (default 80%), can be helpful in cases of false positives
 * 	- complete commenting
 * 	- add "developper" mode
 * 		Disable batch mode and show intermediate steps of analysis in single file mode
 * 	
 * 	Version 2.0
 * 	- add possibility to take fascicle curvature into account
 * 		this is implemented by detecting fascicle orientation in superficial and deep regions and 
 * 		reconstructing a composite fascicle by simple spline fitting of two fascicle segments or 
 * 		by circle fitting. These methods are based on fitting curves across the insertion points and the 
 * 		intersection points between fascicles segments at mid-distance between aponeuroses.
 * 	- detect fascicle orientation based on single (straight fascicle method) or paired 
 * 		(curved fascicles methods) ROIs. This change is possible because the image is now rotated 
 * 		by an angle equal to lower aponeurosis angle and ROIs can be broader.
 * 	- add panoramic mode
 * 		currently, add the possibility to standardise ROI width by selecting a subregion with a 
 * 		width = FoV width * 'x field of view' factor, chosen by user.
 * 		NB: there is no detection of the original FoV width, a width of 5 cm is currently assumed!
 * 	- various changes to filtering method before measurement of fascicle orientation.
 * 		This is important as it means the results are different from previous versions
 * 	- various bug fixes and code improvements.
 * 	
 * 	Version 1.7.1
 * 	- add possibility to extrapolate aponeuroses from 100% or 50% of detected length
 * 	
 * 	Version 1.7
 * 	- add error handling
 * 	- check depencies
 * 	- fix bug that caused analysis to fail when images are flipped and cropping manual
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Script parameters //
///////////////////////
/*									
#@ String(value = "<html><b> SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS </b></html>", visibility="MESSAGE") title
#@ String(value = "NB: Analysis only works with lower aponeurosis orientated horizontally or downwards to the left", visibility="MESSAGE") text0

#@ Boolean (label="Flip image horizontally", value=false, persist=true, description="Required if image displays proximal side to the right") flip
#@ Boolean (label="Panoramic scan", value=false, persist=true, description="panoramic or regular scan") pano
#@ Double(label="x field of view", value = 2.0, min=1.0, max=2.0, stepSize=0.25, persist=false) xFoV
#@ String (label = "Fascicle geometry", choices= {"Straight", "Curved_spline", "Curved_circle"}, style="radioButtonHorizontal", description="Assume straight or curved fascicles. Circle method according to Muramatsu et al JAP 2002") geometry
#@ String (label = "Type of analysis", choices= {"Current file", "Open file", "Open folder"}, style="radioButtonHorizontal", description="Analyse single image or several images") analysis
#@ String (label = "File extension", choices = {".bmp", ".BMP", ".jpg", ".tif", ".png", ".dcm", ""}, description="Required") extension

#@ String(value = "------------------ Single file path -------------------", visibility="MESSAGE") text1
#@ File (label="Select a file", style="file", required = false) inputFile

#@ String(value = "------------------ Directories paths ------------------", visibility="MESSAGE") text2
#@ File (label = "Input directory", style = "directory", required = false) input
#@ File (label = "Output directory", style = "directory", required = false) output

#@ String(value = "------------------ Image cropping ------------------", visibility="MESSAGE") text3	
#@ String(label = " ", choices= {"Automatic (requires identical scan depth)", "Manual"}, style="radioButtonHorizontal", persist=true) cropping

#@ String(value = "--------------- Aponeuroses ---------------", visibility="MESSAGE") text4	
#@ Integer(label = "Tubeness sigma (default: 10)", value = 10, description="Standard deviation of the Gaussian filter. Proportional to aponeurosis thickness") Tsigma
#@ Boolean(label= "Enhance aponeuroses filter", value=false, persist=true, description="run 'Enhance Local Contrast' plugin (CLAHE). Not recommended unless aponeuroses lack contrast") clahe_ap
#@ Integer(label="Length (% of FoV width)", value=80, min=50, max=95, stepSize=5, persist=true, description="Expected length of detected aponeuroses relative to FoV width. default: 80%") apLength
#@ String(label = "Extrapolate from (% of detected aponeurosis length)", choices = {"100%", "50%"}, persist=true, description="") extrapolate_from

#@ String(value = "---------------- Fascicles ----------------", visibility="MESSAGE") text5
#@ Integer(label = "ROI height (% of thickness or 1/2 thickness)", value = 50, min=40, max=90, stepSize=5, persist=true, description="default: 50%") ROIheight
#@ Integer(label = "ROI width (% of ROI width)", value = 60, min=40, max=90, stepSize=10, persist=true, description="default: 60%") ROIwidth
#@ Boolean(label= "Enhance fascicle filter", value=false, persist=true, description="run 'Enhance Local Contrast' plugin (CLAHE). Not recommended unless fascicles lack contrast") clahe_fasc
#@ String(label = "Laplacian of Gaussian (sigma)", choices = {"0", "1", "2", "3", "4", "5", "6", "7", "Test"}, value=1, persist=true, description="Proportional to fascicle thickness: after running the test, choosing the value yielding the greatest angle is recommended") Osigma

#@ String(value = "---------------- Pixel scaling ----------------", visibility="MESSAGE") text6
#@ String(label = "Scaling", choices = {"None", "Automatic (requires metadata)", "Manual"}, value=None, persist=true) scaling
#@ String(label = "Scan depth (cm)", choices = {"3", "4", "5", "6", "7"}, value=5, persist=true) depth

#@ String(value = "-------------------- Other --------------------", visibility="MESSAGE") text7
#@ Boolean(label="Print analysis parameters", value=true, persist=true, description="Add analysis parameters to the results table") param
#@ Boolean(label="Developper mode", value=false, persist=true, description="Disable batch mode and show intermediate steps of analysis") dev
*/

// GLOBAL VARIABLES ------------------------------------------------------------
var is_stack = false;
var image_titles = newArray();
var ALPHA = newArray();
var Lower_aponeurosis_orientation = newArray();
var Pennation_angle = newArray();
var Fascicle_length = newArray();
var Curvature = newArray();
var Thickness = newArray();
var Analysis_duration = newArray();
var Image_cropping = newArray();
var Tubeness_sigma = newArray();
var ROI_h = newArray();
var Thresholding = newArray();
var OrientationJ_sigma = newArray();
var shift_aponeuroses = 0;
var Error = newArray();

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

macro "SMA - Simple Muscle Architecture" {
// description
	plugins = newArray("OrientationJ Measure", "Non-local Means Denoising", "Canny Edge Detector");
	check_dependencies(plugins);
	clear_results_and_roimanager();

	if (analysis == "Current file" || analysis == "Open file") { 
		if (analysis == "Open file") {
			input = File.getDirectory(inputFile);
			file = File.getName(inputFile);
			open_file(input+File.separator+file);
		}
		handle_stacks();
		run("Select None"); 
		run("Remove Overlay");
		scaleFactor = scale(scaling, "");
		if (cropping == "Manual") {
			manual_roi = manual_cropping("");
		} 		
		if (dev == false) {
			setBatchMode(true);
		}		
		singleImageAnalysis();
		
	} else { // process folder
		count = 0;
		countFiles(input);
		n = 0;
		n2 = 1;  // count of processed scans
		scaleFactor = scale(scaling, input);
		if (cropping == "Manual") {
			manual_roi = manual_cropping(input);
		}
		setBatchMode(true);
		processFolder(input);		
		selectWindow("Processed"); 
		run("Close");
	}
	
	
	function processFolder(input) {
		list = getFileList(input);
		list = Array.sort(list);
		run("Text Window...", "name=Processed");
		for (i = 0; i < list.length; i++) {
			print("[Processed]", list[i] +" (scan #"+n2++ +" out of " +count+ ")" +"\n");			
			if (File.isDirectory(input + File.separator + list[i]))
				processFolder(input + File.separator + list[i]);
			if (endsWith(list[i], extension))
				processFile(input, output, list[i]);	
		}
	}
	
	function processFile(input, output, file) {
		open_file(input+File.separator+file);	
		singleImageAnalysis();
	}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// start of singleImageAnalysis function //
///////////////////////////////////////////

	function singleImageAnalysis() {
		
		if (flip == true)
			run("Flip Horizontally", "slice");

		if (pano == true) { // ONLY tested with Philips HD11 and Hologic Mach30
			Wtemp = getWidth; 
			Htemp = getHeight;
			dcmScale = 1 / (metadata_scale() * 10); // Retrieve scaling information in DICOM
			FoV_width = 5;  // 5 cm, arbitrary parameter!!
			pxl_to_FoV = round(5/dcmScale) * xFoV; // x FoV
			makeRectangle(Wtemp-pxl_to_FoV, 0, pxl_to_FoV, Htemp); 		
			run("Crop");
		}
		
		start = getTime(); // Use this to time how long the whole code takes for one image
		title = getTitle();
		IDraw = getImageID();
		W = getWidth; 
		H = getHeight;
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");		
	
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Image Cropping //
		////////////////////		
		if (cropping == "Manual") {
			x = manual_roi[0]; y = manual_roi[1]; width = manual_roi[2]; height = manual_roi[3];
			Le = x+W*0.005;
			Up = y+H*0.01; // shifts cropped window down by 1.5% of its height
			makeRectangle(Le, Up, width*0.99, height*0.98);
		} else {
			run("Duplicate...", " ");
			IDcopy1  = getImageID();
			run("8-bit");
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			run("Convolve...", "text1=[-1 -1 -1 -1 -1\n-1 -1 -1 -1 -1\n-1 -1 24 -1 -1\n-1 -1 -1 -1 -1\n-1 -1 -1 -1 -1\n] normalize");
			run("Median...", "radius=2");
			run("Auto Local Threshold", "method=Median radius=15 parameter_1=0 parameter_2=0 white");
			run("Options...", "iterations=2 count=1 black do=Close");
			run("Analyze Particles...", "size=10000-Infinity add");
			roiManager("Select", 0);
			getSelectionBounds(x, y, width, height);
			roiManager("delete");
			selectImage(IDcopy1); close();
			run( "Select None" );
			selectImage(IDraw);
			Le = x+width*0.005;
			Up = y+height*0.01; // shifts cropped window down by 1% of its height
			makeRectangle(Le, Up, width*0.98, height*0.97);			
		}
	
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Filter and detect aponeuroses //
		///////////////////////////////////

		run("Duplicate...", " ");				
		IDFoV = getImageID();		
		WFoV = getWidth; //Dimensions of the field of view
		HFoV = getHeight;
		setColor(0);
		drawLine(0, 0, WFoV, 0);		
		minLength = apLength/100*WFoV;
		
		// Preprocessing - enhance structures
		run("8-bit");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		run("Subtract Background...", "rolling=50 sliding");	
		run("Non-local Means Denoising", "sigma=15 smoothing_factor=1 auto");
		run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 saturate");
		if (clahe_ap == true) {
			run("Enhance Local Contrast (CLAHE)", "blocksize=36 histogram=256 maximum=4 mask=*None* fast_(less_accurate)");  // necessary for faint aponeuroses
		}
		run("Tubeness", "sigma=Tsigma use"); // currently recommended: 10
		IDvessel1  = getImageID();
		selectImage(IDvessel1);
		run("8-bit");

		// Preprocessing - atenuate fascicles
		// filter out particles with 5 < angle > 50 in quadrant 2
		run("Duplicate...", " ");	
		run("Auto Threshold", "method=Default white");
		run("Median...", "radius=1");					
		bin  = getImageID();
		run("Set Measurements...", "fit redirect=None decimal=3");
		run("Analyze Particles...", "show=Nothing clear record");	
		for (f=0; f<nResults; f++) {
			x = getResult("XStart", f);
		    y = getResult("YStart", f);
		    angle = getResult("Angle", f);
		    area = getResult("Area", f);
			if (angle>130 && angle<175) { 
		        doWand(x,y);
		        run("Clear");
		    }
		}
		run("Select None");
		run("Options...", "iterations=1 count=1 black do=Open");
		run("Create Selection");
		selectImage(bin); close();
		selectImage(IDFoV);							
		run("Restore Selection");
		run("Make Inverse");
		setColor(0);		
		fill();	
		run("Select None");

		// Preprocessing - Threshold from FFT
		run("FFT");
		IDFFT1 = getImageID();
		selectImage(IDvessel1); close();
		selectImage(IDFFT1);
		Wfft = getWidth;
		Hfft = getHeight;
		setThreshold(find_lower_thresh(99.9), 255);
		setOption("BlackBackground", true);
		run("Convert to Mask");		
		setColor(0);
		makeRectangle(Wfft/2+11, 0, Wfft/2, Hfft/2);		
		fill();
		makeRectangle(0, Hfft/2+1, Wfft/2-10, Hfft/2);
		fill();
		run("Inverse FFT");
		IDinvFFT1 = getImageID();
		selectImage(IDFFT1); close();
		selectImage(IDinvFFT1);
		run("Canny Edge Detector", "gaussian=2 low=1 high=7.5"); 			
		run("Analyze Particles...", "size=0-minLength show=Masks");
		IDmask = getImageID(); // Mask of shorter lines
		imageCalculator("Subtract create", IDinvFFT1,IDmask);	
		IDfilteredAp = getImageID();	
		selectImage(IDFoV); close();
		selectImage(IDinvFFT1); close();
		selectImage(IDmask); close();

		// Start detecting position of aponeuroses
		selectImage(IDfilteredAp);
		U = round(0.3*HFoV);
		V = U + 1;
		Upper = lineFinder(0, U, 1);	// The 1 is just a dummy variable so the function can distinguish Upper/Lower
		up_st = newArray(2);
		up_st[0] = Upper[Upper.length-2];
		up_st[1] = Upper[Upper.length-1];
		Upper = Array.slice(Upper, 0, Upper.length-2); // Last 2 elements are line indices
		
		Lower = lineFinder(V, HFoV, 0);
		lo_st = newArray(2);
		lo_st[0] = Lower[Lower.length-2];
		lo_st[1] = Lower[Lower.length-1];
		Lower = Array.slice(Lower, 0, Lower.length-2);
		
		ok = 1;
		selectImage(IDfilteredAp); close();
		
		// warn if aponeurosis was not detected
		err = "";
		if (Upper.length < 10){
			ok = 0;
			if (analysis == "Current file" || analysis == "Open file") {
				exit("Could not detect upper aponeurosis. Please try again with different value of tubeness sigma or try manual cropping");
			} else {
				err = "Could not detect upper aponeurosis. Please try again with different value of tubeness sigma or try manual cropping";
				Error = Array.concat(Error, title +":"+err);
				continue;
			}
		} else if (Lower.length < 10) {
			ok = 0;
			if (analysis == "Current file" || analysis == "Open file") {
				exit("Could not detect lower aponeurosis. Please try again with different value of tubeness sigma or try manual cropping");
			} else {
				err = "Could not detect lower aponeurosis. Please try again with different value of tubeness sigma or try manual cropping";
				Error = Array.concat(Error, title +":"+err);
				continue;
			}
		}
		
		alpha = 0;
		beta = 0;
		theta = 0;
		Lf = 0;
		Th = 0;
		
		if (ok == 1) {
			L = round(0.05*WFoV);
			R = round(0.95*WFoV);
			W2 = R - L;
				
			if (Lower.length < Upper.length) {
				maxL = Upper.length;
				Lower2 = newArray(Upper.length);
				for (t=0; t<Lower.length; t++)
					Lower2[t] = Lower[t];
				i1 = Lower[0]; 
				i2 = Lower[Lower.length-1]; 
				islope = (i2 - i1) / Lower.length;
				for (p=Lower.length; p<Upper.length; p++)
					Lower2[p] = Lower2[p-1] + islope;
				Lower = Array.copy(Lower2);
			} else if (Upper.length < Lower.length) {
				maxL = Lower.length;
				Upper2 = newArray(Lower.length);
				for (t=0; t<Upper.length; t++)
					Upper2[t] = Upper[t];
				i1 = Upper[0]; 
				i2 = Upper[Upper.length-1]; 
				islope = (i2 - i1) / Upper.length;
				for (p=Upper.length; p<Lower.length; p++)
					Upper2[p] = Upper2[p-1] + islope;
				Upper = Array.copy(Upper2);
			} else
				maxL = Upper.length; // If both aponeuroses are the same length, we can use either length
			
			x_upp = newArray(Upper.length); // Create x values for overlay
			for(t=0; t<Upper.length; t++){
				Upper[t] = Upper[t] + Up;
				x_upp[t] = up_st[0] + t + Le;
			}
			x_low = newArray(Lower.length);
			for(t=0; t<Lower.length; t++){
				Lower[t] = Lower[t] + Up;
				x_low[t] = lo_st[0] + t + Le;
			}
			
			// Plot the curves on top of the existing image	
			selectImage(IDraw);
			applyOverlay(Upper, x_upp);
			applyOverlay(Lower, x_low);

			// select part of aponeuroses to extrapolate from
			upper_idx1 = 0;
			lower_idx2 = x_low.length-1;
			if (extrapolate_from == "100%") {
				upper_idx2 = x_upp.length-1;
				lower_idx1 = 0;
			} else {  // 50%
				mid_upper = middle_coordinates(x_upp, Upper);
				mid_lower = middle_coordinates(x_low, Lower);
				upper_idx2 = mid_upper[2];
				lower_idx1 = mid_lower[2];
			}
			
			//Upper aponeurosis line (y_1 = b_1 * x + a_1)
			Fit.doFit("Straight Line", Array.slice(x_upp, upper_idx1, upper_idx2), Array.slice(Upper, upper_idx1, upper_idx2));
			a_1 = Fit.p(0); b_1 = Fit.p(1);
			betaUp = atan(b_1)* (180/PI); //angle upper aponeurosis
			//Lower aponeurosis line (y_2 = b_2 * x + a_2)
			Fit.doFit("Straight Line", Array.slice(x_low, lower_idx1, lower_idx2), Array.slice(Lower, lower_idx1, lower_idx2));
			a_2 = Fit.p(0); b_2 = Fit.p(1);
			betaDeep = atan(b_2)* (180/PI); //angle deep aponeurosis

			// rename aponeuroses coordinates arrays 
			UAx = Array.copy(x_upp); //upper aponeurosis x values
			UAy = Array.copy(Upper); //upper aponeurosis y values
			LAx = Array.copy(x_low); //lower aponeurosis x values
			LAy = Array.copy(Lower); //lower aponeurosis y values

			// calculate mid-distance line between aponeuroses and get equation (2nd order polynom)
			MAx = Array.copy(UAx);
			MAy = newArray(UAy.length);
			for (m=0; m<UAx.length; m++){
				MAy[m]=UAy[m]+(LAy[m]-UAy[m])*0.5; // @ 0.5 * thickness from top aponeurosis
			}
			//Mid-distance line (y_m = b_m * x + a_m)      
			b_m = (MAy[0] - MAy[MAx.length-1]) / (MAx[0] - MAx[MAx.length-1]); //slope based on whole aponeurosis
			a_m = MAy[0] - b_m * MAx[0];
			betaM = atan(b_m)* (180/PI); //angle mid line
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Measure fascicle orientation	//
		//////////////////////////////////	
			
		// Select ROI between aponeuroses
		// retrieve statistics from aponeuroses
		Array.getStatistics(UAy, minUAy, maxUAy, meanUAy, stdUAy); 
		Array.getStatistics(LAy, minLAy, maxLAy, meanLAy, stdLAy); 			
		mindist = minLAy-maxUAy; //smallest distance between aponeuroses

		// store aponeuroses arrays as ROIs
		makeSelection("polyline", UAx, UAy);
		roiManager("add");
		makeSelection("polyline", LAx, LAy);
		roiManager("add");

		// rotate image by angle equal to lower aponeurosis orientation
		run("Duplicate...", " ");
		IDrawR  = getImageID();
		betaDeepR = -betaDeep;
        run("Rotate... ", "angle=betaDeepR grid=1 interpolation=Bicubic");		

		// rotate stored aponeuroses ROIs
		roiManager("Select", 0);
		run("Rotate...", "  angle=betaDeepR");
		roiManager("update");
		roiManager("Select", 1);
		run("Rotate...", "  angle=betaDeepR");
		roiManager("update");

		// retrieve rotated aponeuroses coordinates
		roiManager("Select", 0);
		Roi.getCoordinates(UAxR, UAyR);
		roiManager("Select", 1);
		Roi.getCoordinates(LAxR, LAyR);

		// statistics on rotated aponeuroses
		Array.getStatistics(UAyR, minUAyR, maxUAyR, meanUAyR, stdUAyR); 
		Array.getStatistics(LAyR, minLAyR, maxLAyR, meanLAyR, stdLAyR); 
		mindistR = minLAyR-maxUAyR;

		// make ROI(s) for detection of fascicle angle
		if (geometry == "Curved_spline" || geometry == "Curved_circle") {
			makeRectangle(UAxR[0], maxUAyR+5, (UAxR.length-1)*ROIwidth/100, mindistR*ROIheight*0.01); //superficial			
			roiManager("add");
		}
		if (geometry == "Straight") {
			deep_roi_x = UAxR[(UAxR.length-1)*0.5] - (UAxR.length-1)*ROIwidth/100*0.5;
		} else {
			deep_roi_x = UAxR[(UAxR.length-1)*(1-ROIwidth/100)];
		}
		makeRectangle(deep_roi_x, minLAyR-mindistR*(ROIheight+5)*0.01, (UAxR.length-1)*ROIwidth/100, mindistR*ROIheight*0.01); //deep			
		roiManager("add");				
		roiManager("show all");

		// Select fascicle area ROI(s)
		n = roiManager("count"); 			
		if (n > 0) {
			roiManager("Show All");
			cont = 1; 
		}
		if (cont > 0) {
			store = newArray(n);										
			for (i=2; i<n; i++)  {
				roiManager("select", i);						
				run("Duplicate...", " ");
				IDROI  = getImageID();
				
				// Preprocessing - enhance structures
				run("32-bit");
				run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
				run("Subtract Background...", "rolling=50 sliding");
				run("Non-local Means Denoising", "sigma=15 smoothing_factor=1 auto");
				if (clahe_fasc == true) {
					run("Enhance Local Contrast (CLAHE)", "blocksize=54 histogram=256 maximum=32 mask=*None* fast_(less_accurate)");
				}				
				run("Median...", "radius=2");	
				run("Tubeness", "sigma=2 use");
				run("8-bit");
				
				// Preprocessing - atenuate fascicles
				// filter out particles with angle <50 or > 5 in quadrant 2, or with an area < 50 pxls
				run("Auto Threshold", "method=Default white");
				run("Median...", "radius=1");
				bin  = getImageID();
				run("Set Measurements...", "fit redirect=None decimal=3");
				run("Analyze Particles...", "show=Nothing clear record");
				for (f=0; f<nResults; f++) {
				    x = getResult("XStart", f);
				    y = getResult("YStart", f);
				     angle = getResult("Angle", f);
				    area = getResult("Area", f);
				    if (angle<130 || angle>175 || area<50) {
				        doWand(x,y);
				        run("Clear");
				    }
				}
				run("Select None");
				run("Options...", "iterations=1 count=1 black do=Open");
				run("Create Selection");
				selectImage(bin);	
				close ();
				selectImage(IDROI);							
				run("Restore Selection");
				run("Make Inverse");
				setColor(0);		
				fill();	
				run("Select None");		
				run("FFT");
				IDFFT2 = getImageID();
				selectImage(IDFFT2);									
				setThreshold(find_lower_thresh(99.7), 255);  // used to be 99.82
				run("Convert to Mask");
				run("Inverse FFT");
				IDinvFFT2 = getImageID();
				selectImage(IDFFT2);
				close ();
				selectImage(IDinvFFT2);
				run("8-bit");
				n_frangi = 3;  // number of times to run the Frangi filter
				for (f = 0; f < n_frangi; f++) {
					temp = getImageID();
					selectImage(temp);
					run("Frangi Vesselness (imglib, experimental)", "number=1 minimum=1.000000 maximum=1.000000");
					selectImage(temp); close();
				}
				IDvessel2  = getImageID();
				selectImage(IDvessel2);	

				// run OrientationJ plugin
				run("Select All");
				if (Osigma == "test") {
					for (o = 0; o < 8; o++) {
						run("OrientationJ Measure", "sigma="+o);
					}
					exit();
				} else {
					run("OrientationJ Measure", "sigma="+Osigma);
				}
				selectImage(IDvessel2);
				close ();
				
				// retrieve angle values from OrientationJ log window
				table = getInfo("log"); 					
				lines = split(table, "\n"); 
				headings = split(lines[0], "\t"); 
				values = split(lines[1], "\t");
				selectWindow("Log");
				run("Close");
				store[i] = values[6];
				store[i] = toString(store[i]);
				store[i] = replace(store[i], ",", ".");
				store[i] = parseFloat(store[i]);

				selectImage(IDROI);
				close ();	
			}
		}
		selectImage(IDrawR);
		close ();
		roiManager("reset");			

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Calculate fascicle length, pennation angle and thickness //
		//////////////////////////////////////////////////////////////
			
		if (geometry == "Curved_spline" || geometry == "Curved_circle") {
		
			alphaDeep = -store[3]; // with angle due to rotation
			alphaUp = -store[2];		
		    
			//Pennation angle		
			thetaDeep = alphaDeep + betaDeep; //orientation of lower part of fasc. relative to deep apon.
			thetaUp = alphaUp + betaDeep; //orientation of upper part of fasc. relative to deep apon.		

			//Composite fascicle lines (y_3 = b_3 * x + a_3)
			b_3d = tan(thetaDeep * (PI / 180)); // slope for the deeper region
			b_3u = tan(thetaUp * (PI / 180)); // slope for the upper region			
			a_3d = LAy[LAy.length-1] - b_3d * LAx[LAx.length-1]; 
			a_3u = MAy[MAy.length-1] - b_3u * MAx[MAx.length-1];		

			//Fascicle parts intersection coordinates	
			IxFasc = MAx[MAx.length/2];	// intersection placed on point in the middle of mid line
			IyFasc = MAy[MAy.length/2];			

			//Upper fascicle part - upper aponeurosis intersection coordinates
			tempxUp = (a_1 - a_3u) / (b_3u - b_1);
			tempyUp = a_1 + b_1 * tempxUp;			
			
			//Offset upper fascicle part  (y_3 = b_3 * x + a_3)
			b_Uo = tan(thetaUp * (PI / 180));						
			a_Uo = IyFasc - b_Uo * IxFasc; 
				
			//Offset coordinates on upper aponeurosis for upper part of fascicle 
			Oux = (a_1 - a_Uo) / (b_Uo - b_1);
			Ouy = a_1 + b_1 * Oux;	
		
			//Deeper fascicle part - lower aponeurosis intersection coordinates
			tempxDeep = (a_2 - a_3d) / (b_3d - b_2);
			tempyDeep = a_2 + b_2 * tempxDeep;	
			
			//Offset lower fascicle part  (y_3 = b_3 * x + a_3)
			b_Do = tan(thetaDeep * (PI / 180));							
			a_Do = IyFasc - b_Do * IxFasc; 	
			
			//Offset coordinates on lower aponeurosis for lower part of fascicle 
			Olx = (a_2 - a_Do) / (b_Do - b_2);
			Oly = a_2 + b_2 * Olx;										

			// From https://forum.image.sc/t/how-do-you-measure-the-angle-of-curvature-of-a-claw/2056/16			
			// get center of circle passing through fascicles intersection points and mid-depth fascicle point
			x_centre = (b_Uo*b_Do*(-Ouy+Oly)+b_Uo*(IxFasc+Olx)-b_Do*(IxFasc+Oux))/(2*(b_Uo-b_Do));
			y_centre = -1/b_Uo*(x_centre-(IxFasc+Oux)/2)+(IyFasc+Ouy)/2;
			
			// angle, slope and intercept of segments between circle centre and fascicle intersections, and angle (gamma) between them
			b_centreUp = (y_centre - Ouy) / (x_centre - Oux);
			a_centreUp = -(y_centre - Ouy) / (x_centre - Oux) * x_centre + y_centre;	
			ang_Up = atan(b_centreUp)* (180/PI); //angle superficial segment				
			b_centreDeep = (y_centre - Oly) / (x_centre - Olx);
			a_centreDeep = -(y_centre - Oly) / (x_centre - Olx) * x_centre + y_centre;
			ang_Deep = atan(b_centreDeep)* (180/PI); //angle deep segment	
			
			tan_gamma = abs((b_centreUp-b_centreDeep)/(1+b_centreDeep*b_centreUp));
			gamma = atan(tan_gamma)* (180/PI);
			
			//radius, axes, curvature and fascicle length
			r = sqrt((x_centre-Oux)*(x_centre-Oux) + (y_centre-Ouy)*(y_centre-Ouy)); //radius (to upper insertion)
			r2 = sqrt((x_centre-Olx)*(x_centre-Olx) + (y_centre-Oly)*(y_centre-Oly)); //radius (to lower insertion)				
			if (y_centre > 0) { //fascicle curved towards upper aponeurosis
				Crv = 1/r; //curvature
			} else {
				Crv = -1/r;
			}
			
			if (geometry == "Curved_circle") { 
				Lf = 2*PI*r*(gamma/360);
				c_x = newArray(round(gamma+1));
				c_y = newArray(round(gamma+1));
				if (y_centre > 0) {
					for (i = 0; i < round(gamma)+1; i++) {
					c_x[i]  =  x_centre + r * cos((abs(ang_Deep)+i)*PI/180);
					c_y[i]  =  y_centre - r * sin((abs(ang_Deep)+i)*PI/180); // "-" because the Y origin points downwards
					}
				} else {
					for (i = 0; i < round(gamma)+1; i++) {
					c_x[i]  =  x_centre - r * cos((abs(ang_Up)+i)*PI/180);
					c_y[i]  =  y_centre + r * sin((abs(ang_Up)+i)*PI/180); // "-" because the Y origin points downwards
					}
				}
				makeSelection(7, c_x, c_y); // make circle arc
				roiManager("add");
														
			} else { // Spline
				makeSelection("polyline", newArray(Oux,IxFasc,Olx), newArray(Ouy,IyFasc,Oly));
				run("Fit Spline");
				run("Interpolate");
				List.setMeasurements;
				Lf = List.getValue("Length");				
				roiManager("add");
			} 
		
		} else { //if fascicles are analysed as straight lines
			alphaDeep = -store[2];			
		    
			//Pennation angle
			theta = alphaDeep + betaDeep;	

			//Composite fascicle line  (y_3 = b_3 * x + a_3)
			b_3 = tan(theta * (PI / 180));				
			a_3 = LAy[LAy.length-1] - b_3 * LAx[LAx.length-1]; //Assumption: the upmost right point on the lower aponeurosis
				
			//Fascicle-upper aponeurosis intersection coordinates
			Ix = (a_1 - a_3) / (b_3 - b_1);
			Iy = a_1 + b_1 * Ix;

			//Offset coordinates on lower aponeurosis
			Hx1 = (Ix + b_2 * Iy - b_2 * a_2) / ((b_2 * b_2) + 1);
			Hy1 = b_2 * ((Ix + b_2 * Iy - b_2 * a_2) / ((b_2 * b_2) + 1)) + a_2;	
			Olx = LAx[LAx.length-1] + ((LAx[0]-Hx1)/2);	
			Oly = LAy[LAy.length-1] + ((LAy[0]-Hy1)/2);			

			//Offset composite fascicle line  (y_3 = b_3 * x + a_3)
			b_o = tan(theta * (PI / 180));
			a_o = Oly - b_o * Olx; 
			
			//Offset coordinates on upper aponeurosis
			Oux = (a_1 - a_o) / (b_o - b_1);
			Ouy = a_1 + b_1 * Oux;	
			//Offset coordinates on lower aponeurosis
			Olx = (a_2 - a_o) / (b_o - b_2);
			Oly = a_2 + b_2 * Olx;	

			//unscaled fascicle length
			Lf = sqrt((Olx-Oux)*(Olx-Oux) + (Oly-Ouy)*(Oly-Ouy)); //fascicle length
			Crv = NaN; //no curvature
			r = NaN;
	
			makeLine(Olx, Oly, Oux, Ouy);
			roiManager("Add");	
		}

		//adjust fascicle and ROIs display
		Omx = (Olx+Oux)/2; //x coordinate of mid fascicle
		mFoVx = UAx[(UAx.length)/2]; //x coordinate of mid FoV			
		if (Oux<0 || Olx>W) {
			n = roiManager("count");
			roiManager("select", n-1);
			roiManager("translate", -(Omx-mFoVx), 0);
			roiManager("Update");
			if (geometry == "Straight") { 
				getLine(Nx1, Ny1, Nx2, Ny2, NlineWidth); // find start and end coordinates of straight line
			} else {					
				Roi.getCoordinates(Spl_xpoints, Spl_ypoints) // find start and end coordinates of spline
				Nx2 = Spl_xpoints[0];	
				Ny2 = Spl_ypoints[0];	
				Nx1 = Spl_xpoints[Spl_xpoints.length-1];	
				Ny1 = Spl_ypoints[Spl_ypoints.length-1];			
			}
			left = abs(Nx2);
			right = Nx1-W;			
			if (Nx2<0 && Nx1<W) {
				Wnew=W+left; 
				run("Canvas Size...", "width=Wnew height=H position=Center-Right zero"); 
			} else if (Nx2>0 && Nx1>W) {
				Wnew=W+right; 
				run("Canvas Size...", "width=Wnew height=H position=Center-Left zero");          			
			} else if (Nx2<0 && Nx1>W && left>right) {
				shift_aponeuroses = left;
				Wnew=W+left*2;
				run("Canvas Size...", "width=Wnew height=H position=Center zero"); 
				for (i=0; i<n; i++)  {
					roiManager("select", i);
				    roiManager("translate", left, 0);
					roiManager("Update");
				}  
			} else if (Nx2<0 && Nx1>W && left<right) {
				shift_aponeuroses = right;
				Wnew=W+right*2;
				run("Canvas Size...", "width=Wnew height=H position=Center zero"); 	 
				for (i=0; i<n; i++)  {
					roiManager("select", i);
					roiManager("translate", right, 0);
					roiManager("Update");
				}	       			
			}
		}
			
		//scaled fascicle length and curvature
		if (scaling == "Automatic (requires metadata)" || scaling == "Manual") {
			Lf = Lf / scaleFactor;
			scaleFactor2 = scaleFactor/1000;
			Crv = 1/(r/scaleFactor2); //different scaling factor to get the curvature in m-1
		}
		//Thickness (temporary). 
		sta = minOf(x_upp.length, x_low.length);
		Th = newArray(sta);
		sum = 0;
		for (i=0; i<Th.length; i++){
			Th[i] = Lower[i] - Upper[i];
			sum += Th[i];
		}
		
		if (scaling == "Automatic (requires metadata)" || scaling == "Manual") {
			Th = (sum/Th.length)/scaleFactor;	
		} else {
			Th = sum/Th.length;
		}

		// Analysis duration
		Duration = (getTime() - start) / 1000;
			
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Save variables and figures //
		////////////////////////////////
		
		// Save variables of interest
		if (geometry == "Curved_spline" || geometry == "Curved_circle") {		
			image_titles = Array.concat(image_titles, title);
			Lower_aponeurosis_orientation = Array.concat(Lower_aponeurosis_orientation, betaDeep);
			Pennation_angle = Array.concat(Pennation_angle, alphaDeep);
			Fascicle_length = Array.concat(Fascicle_length, Lf);
			Curvature = Array.concat(Curvature, Crv);
			Thickness = Array.concat(Thickness, Th);
			Analysis_duration = Array.concat(Analysis_duration, Duration);
			Image_cropping = Array.concat(Image_cropping, cropping);
			Tubeness_sigma = Array.concat(Tubeness_sigma, Tsigma);
			ROI_h = Array.concat(ROI_h, ROIheight);
			OrientationJ_sigma = Array.concat(OrientationJ_sigma, Osigma);
		} else {  // process folder
			image_titles = Array.concat(image_titles, title);
			Lower_aponeurosis_orientation = Array.concat(Lower_aponeurosis_orientation, betaDeep);
			Pennation_angle = Array.concat(Pennation_angle, alphaDeep);		
			Fascicle_length = Array.concat(Fascicle_length, Lf);
			Thickness = Array.concat(Thickness, Th);
			Analysis_duration = Array.concat(Analysis_duration, Duration);
			Image_cropping = Array.concat(Image_cropping, cropping);
			Tubeness_sigma = Array.concat(Tubeness_sigma, Tsigma);
			ROI_h = Array.concat(ROI_h, ROIheight);
			OrientationJ_sigma = Array.concat(OrientationJ_sigma, Osigma);
		} 
	    
	    // save figures with overlays
	    selectImage(IDraw);
	    roiManager("Set Color", "yellow");
		roiManager("Set Line Width", 4);
		run("From ROI Manager");
	    if (analysis == "Open folder") {
			run("Flatten", "stack");
			saveAs("tiff", output + File.separator + file);
			roiManager("reset");
			close("*");
		} else {
			selectImage(IDraw);
			run("Flatten", "slice");
			selectImage(IDraw);
		}
		
		// Convert 0 values to NaN for final display
		ALPHA = arrTidy(ALPHA);
		Lower_aponeurosis_orientation = arrTidy(Lower_aponeurosis_orientation);
		Pennation_angle = arrTidy(Pennation_angle);
		Fascicle_length = arrTidy(Fascicle_length);
		Curvature = arrTidy(Curvature);
		Thickness = arrTidy(Thickness);		
	
	} //end of singleImageAnlysis function
	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Display results //
	/////////////////////
	function display_results() { 
	// function description
		if (param == false) {
			Array.show("Results",
			image_titles,
			Lower_aponeurosis_orientation,
			Pennation_angle,
			Fascicle_length,
			Curvature,
			Thickness,
			Analysis_duration);
		} else {
			Array.show("Results",
			image_titles,
			Lower_aponeurosis_orientation,	
			Pennation_angle,
			Fascicle_length,
			Curvature,
			Thickness,
			Analysis_duration,
			Image_cropping,
			Tubeness_sigma,
			ROI_h,
			Thresholding,
			OrientationJ_sigma);
		}
	}
	display_results();
	
	function displayStack(output) {
		stack1=output+File.separator; 
		run("Image Sequence...", "open=&stack1 sort use");
	}

	if (analysis == "Open folder") {  
		list = getFileList(output);
		displayStack(output);
    }
	if (Error.length > 0) {
		Array.show("Images not processed", Error);
	}
	if (dev == false) {
		setBatchMode(false);	
	}
} //end of macro	 
	  

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Help functions //
////////////////////

function check_dependencies(deps) { //TODO
// function description
	missing = newArray();
	List.setCommands; 
	for (i=0; i<deps.length; i++) {
		if (List.get(deps[i])=="") { 
			missing = Array.concat(missing, deps[i]);
		}
	}
	if (missing.length > 0) {
		exit("SMA requires the following plugins:"
			+"\n "
			+"\n- OrientationJ,"
			+"\n- Non-local Means Denoising,"
			+"\n- Canny Edge Detector"
			+"\n "
			+"\nPlease chack that the update sites 'BIG-EPFL', 'Biomedgroup' and 'BioVoxxel' are selected."
			+"\n\"https://imagej.net/Following_an_update_site.html\"");
	}
}


function clear_results_and_roimanager() { 
// function description
	if (isOpen("Results")) { //Close results window (obtained from previous analysis, should be kept out of loop for stacks)
		selectWindow("Results");
	    run("Close");
	}
	if (isOpen("ROI Manager")) {
		selectWindow("ROI Manager");
	 	run("Close");
	}
}


function countFiles(input) {
  list = getFileList(input);
  for (i=0; i<list.length; i++) {
      if (endsWith(list[i], "/"))
          countFiles(""+input+list[i]);
      else
          count++;
  }
}


function open_file(path)	{
// Open single DICOM file
//	if (extension == "" || extension == ".dcm") {
	splitted = split(path, ".");
	if (path.endsWith(".dcm") || splitted.length == 1) {
		print(File.getName(path));
		run("Bio-Formats","open=path autoscale color_mode=Default open_files view=Hyperstack stack_order=XYCZT");
		unstack_channels();
	} else if (path.endsWith(".avi")) {
		run("Movie (FFMPEG)...", "choose=[path] first_frame=0 last_frame=-1");
	} else {
		open(path);
	}
}


function unstack_channels() { 
// convert image stack to RGB image
	title = getTitle();
	getDimensions(im_width, im_height, im_channels, im_slices, im_frames);
	if (im_channels > 1) {
		run("Stack to RGB");
		run("8-bit");
		selectWindow(title); close();
	}
}


function scale(scaling_mode, folder_path) { 
// function description
	if (folder_path.length > 0) {
		list = getFileList(folder_path);
		open_file(folder_path+File.separator+list[0]);
	}
	if (scaling_mode == "Automatic (requires metadata)") {
		scaleFactor = metadata_scale();
	} else if (scaling_mode == "Manual") {
		scaleFactor = userscale();
	}
	if (folder_path.length > 0) {
		close ();
	}
	return scaleFactor;
}


function manual_cropping(folder_path) {
// function description
	if (folder_path.length > 0) {
		list = getFileList(folder_path);
		open_file(folder_path+File.separator+list[0]);
		if (flip == true) {
			run("Flip Horizontally");
		}
	} else if (is_stack == true) {
		run("Duplicate...", " ");
		if (flip == true) {
			run("Flip Horizontally");
	}
	setTool("rectangle"); 
	waitForUser("Select area. Click OK when done");
	Roi.getBounds(x, y, width, height);
	if (folder_path.length > 0 || is_stack == true) {
		close ();
	}
	return newArray(x, y, width, height);		
}


function getNumericTag(tag) { // get value from DICOM file. Source: https://imagej.nih.gov/ij/macros/ListDicomTags.txt
    value = getTag(tag);
    if (value=="") return NaN;
    index3 = indexOf(value, "\\");
    if (index3>0)
      value = substring(value, 0, index3);
    value = 0 + value; // convert to number
    return value;
  }
// This function returns the value of the specified 
// tag  (e.g., "0010,0010") as a string. Returns "" 
// if the tag is not found.
function getTag(tag) {
  info = getImageInfo();
  index1 = indexOf(info, tag);
  if (index1==-1) return "";
  index1 = indexOf(info, ":", index1);
  if (index1==-1) return "";
  index2 = indexOf(info, "\n", index1);
  value = substring(info, index1+1, index2);
  return value;
}
  
   
function arrTidy(Arr) {
	for (m=0; m<Arr.length; m++) {
		if (Arr[m] == 0)
			Arr[m] = NaN;
	}
	return Arr;
}


function find_lower_thresh(excluded_percent) { 
// return lower value used to threshold a FFT image, based on percentage of pixels that should be excluded
	percentage = excluded_percent; // percentage of FFT pixels to threshold out						
	nBins = 256; 
	resetMinAndMax(); 
	getHistogram(pxlValues, histCounts, nBins); 
	// find culmulative sum 
	nPixels = 0; 
	for (f = 0; f<histCounts.length; f++) 
		nPixels += histCounts[f]; 
	nBelowThreshold = nPixels * percentage / 100; 
	sum = 0; 
	for (f = 0; f<histCounts.length; f++) { 
		sum = sum + histCounts[f]; 
		if (sum >= nBelowThreshold) { 
			lowThresh = pxlValues[f]; // lower threshold corresponding to % of highest GV pixels in FFT   
	    	f = 99999999;//break 
	  	} 
	}
	return lowThresh;
}


function lineFinder(in1, in2, in3) { 
//in1: vert start point of line, in2: vert end point of line, in3: 1 for upper line, 0 for lower
	L = round(0.05*WFoV); R = round(0.95*WFoV); M = 0.5*WFoV;	
	makeLine(M, in1, M, in2); // Search from middle of image
	profile = getProfile();
	hold = newArray(1);
	for (j=0; j<profile.length; j++) // First count how many 'lines' there are
		if (profile[j] > 50)  // 50 is a threshold- pixel intensities above 50 are white
			hold = append(hold, j + in1, 0); // Store the indices of each possible line in hold
	hold = Array.slice(hold, 1, hold.length);
	line = newArray(1);
	ind1 = 0;
	ind2 = 1;
	
	if (in3 == 1)
		a = hold.length - 1;
	else
		a = 0;
	found = 0; 
	while (found != 1 && a >= 0 && a < hold.length) {
		line[0] = hold[a]; // Start with right half of line
		if(in3 == 1)
			makeLine(M+1, line[0]-21, M+1, line[0]-1); // Check just above/below each 'line' for a parallel line
		else
			makeLine(M+1, line[0]+1, M+1, line[0]+20);
		pr = getProfile();
		sum = 0;
		for(c=0; c<pr.length; c++)
			sum += pr[c];
		if(sum < 50)
			sum += 1; // Do nothing (this is not the right line, try the next one)
		else {
		
			b = line[0];
			for(s=M+1; s<R; s++) {
				makeLine(s, b-1, s, b+1); // Assume aponeurosis location varies by no more than +/-1 per pixel
				prof = getProfile();
				
				for(t=0; t<prof.length; t++) {					
					if(prof[t] > 50) {
						b = (b - 1) + t;
						line = append(line, b, 0);  
						t = prof.length; 
						ind2 = s; } }
					}
					
			b = line[0];
			// Find left half of line
			for(s=M; s>L; s--) {  // WAS M-1, might change back
				makeLine(s, b-1, s, b+1); // Assume aponeurosis location varies by no more than +/-1 per pixel
				prof = getProfile();
				
				for(t=0; t<prof.length; t++) {
					if(prof[t] > 50) {
						b = (b - 1) + t; 
						line = append(line, b, 1); 
						t = prof.length; 
						ind1 = s; } } 
								}
	}
	if (line.length > 2) { // If a suitable line is found, stop the loop
		a = hold.length;
		found = 1; }
	else if (in3 == 1 && a-1 >= 0) // Or increment and keep looking for a suitable line
		a -= 1;
	else if (in3 == 0 && a+1 < hold.length)
		a += 1;
	else	// Otherwise stop the loop and move to the next section
		a = hold.length;
	}
	// If no suitable line is found, find the most proximal (deep aponeurosis) or distal line
		if(found == 0) { 
			for(a=0; a<hold.length; a++){
				line[0] = hold[a];
				b = line[0];
				for(s=M+1; s<R; s++) {
				makeLine(s, b-1, s, b+1); // Assume aponeurosis location varies by no more than +/-1 per pixel
				prof = getProfile();
				// Right half of line
				for(t=0; t<prof.length; t++) {					
					if(prof[t] > 50) {
						b = (b - 1) + t;
						line = append(line, b, 0);  
						t = prof.length;
						ind2 = s; } }
					}
				b = line[0];

			// Find left half of line
			for(s=M-1; s>L; s--) {
				makeLine(s, b-1, s, b+1); // Assume aponeurosis location varies by no more than +/-1 per pixel
				prof = getProfile();
				
				for(t=0; t<prof.length; t++) {
					if(prof[t] > 50) {
						b = (b - 1) + t; 
						line = append(line, b, 1); 
						t = prof.length;
						ind1 = s; } } 
								}
			if(line.length > 0.60*WFoV)
				a = hold.length;
		}
		}
	line = append(line, ind1, 0); // Need the indices of the start/end of the accepted line
	line = append(line, ind2, 0);
	return line; 
	}


function compDiff() {
	oz = newArray(p.length-1);
	for(t=0; t<p.length-1; t++){
		oz[t] = p[t+1] - p[t]; // Compute diff pixel by pixel
		if(oz[t] != 0) oz[t] = 1; // Convert to binary
	}
	return oz;
}


function applyOverlay(cInput, x) {	
	run("Select All");
	makeSelection("polyline", x, cInput);
	run("Fit Spline");
	run("Overlay Options...", "stroke=green width=4 set");
	run("Add Selection...");
	run("Select None"); }


function append(arr, value, place) {
// arr = array to append to, value = value to append, place = 1 for start, 0 for end
    arr2 = newArray(arr.length+1);

	if (place == 1) { // Append to start
		arr2[0] = value;
		for (i=0; i<arr.length; i++)
        	arr2[i+1] = arr[i];
     	return arr2; }

    else
    	{
     	for (i=0; i<arr.length; i++) // Append to end
        	arr2[i] = arr[i];
     	arr2[arr.length] = value;
     	return arr2; }
  	}


function index(a1, a2, value, condition) { 
	for (i=0; i<a1.length; i++) 
	  if ((a1[i]>=round(value) || a1[i]==round(value+1) || a1[i]==round(value-1)) && (a2[i]==round(condition) || a2[i]==round(condition+1) || a2[i]==round(condition-1))) return i; 
	return -1; 
} 


function middle_coordinates(xPoints, yPoints) {
// 
	nPoints = xPoints.length - 1;
	lineLength = 0;
	for (i=1; i<nPoints; i++) {
	  xDistance = xPoints[i] - xPoints[i-1];
	  yDistance = yPoints[i] - yPoints[i-1];
	  segmentLength = sqrt(xDistance*xDistance + yDistance*yDistance);
	  lineLength += segmentLength;
	}
	halfLineLength = lineLength / 2;
	runningLength = 0;
	for (i=1; i<nPoints; i++) {
	  xDistance = xPoints[i] - xPoints[i-1];
	  yDistance = yPoints[i] - yPoints[i-1];
	  segmentLength = sqrt(xDistance*xDistance + yDistance*yDistance);
	  runningLength += segmentLength;
	  if (runningLength >= halfLineLength) {
	    midX = xPoints[i-1] + (xPoints[i] - xPoints[i-1]) * (halfLineLength - (runningLength - segmentLength)) / segmentLength;
	    midY = yPoints[i-1] + (yPoints[i] - yPoints[i-1]) * (halfLineLength - (runningLength - segmentLength)) / segmentLength;
	    midSegment = i - 1;
	    break;
	  }
	}
	return newArray(midX, midY, midSegment);
}


function metadata_scale() { 
// function description
	metadata = getMetadata("Info");
	if (metadata.length > 0) {
		dcmScale = parseFloat(getInfo("Physical Delta X"));  // original scale in cm
	} else {
		exit("No metadata data were found and the scaling can't be done automatically");
	}
	return 1/dcmScale/10;  // scaling factor
}


function userscale() { 
// function description
	run("Select None");
	setTool("line");
	waitForUser("Select scaling line. Click OK when done");
	getLine(x1, y1, x2, y2, lineWidth);
	lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
	return lineLength/(parseInt(depth)*10);
}


function handle_stacks() { 
// function description
	if (nSlices == 1)  // not a stack
		return;
	is_stack = true;
	
	Dialog.createNonBlocking("Image stack");
	choices = newArray("Current frame only", "Multiple frames");
	Dialog.addRadioButtonGroup("Process:", choices, 2, 1, "Current frame only");
	Dialog.addSlider("start", 1, nSlices, 1);
	Dialog.addSlider("end", 1, nSlices, nSlices);
	Dialog.addMessage("Alternatively, enter comma-separated frame numbers");
	Dialog.addString("Discrete frames", "", 4);
	Dialog.show();
	target = Dialog.getRadioButton();
	start = Dialog.getNumber();
	end = Dialog.getNumber();
	slices_string = Dialog.getString();

	if (target == "Current frame only") {
		end = 1;
	}
	if (slices_string != "") {
		slices = split(slices_string, ",");
		for (i = 0; i < slices.length; i++) {
			slices[i] = parseInt(slices[i]);
		}
		start = 0; end = slices.length - 1;
	}
	scaleFactor = scale(scaling, "");
	if (cropping == "Manual") {
		manual_roi = manual_cropping("");
	} 
	if (dev == false) {
		setBatchMode(true);
	}
	for (slice_number=start; slice_number<end+1; slice_number++) {
		if (slices_string != "") {
			Stack.setSlice(slices[slice_number]);
		} else {
			Stack.setSlice(slice_number);
		}
		Overlay.hide
		run("Select None"); 
		singleImageAnalysis();
		roiManager("reset");
		display_results();
	}
	Overlay.show
	exit;
}
