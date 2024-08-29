///////////////////////////////////////////////
// SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS //
///////////////////////////////////////////////
/*
*   Reference:
*   Seynnes OR, Cronin NJ. Simple Muscle Architecture Analysis (SMA): An ImageJ macro tool to automate measurements in B-mode ultrasound scans. PLoS One. 2020;15: e0229034. doi:10.1371/journal.pone.0229034
*   Hosted at: 
*   https://github.com/oseynnes/SMA
*   
*	Requires OrientationJ plugin v2.0.3 (http://bigwww.epfl.ch/demo/orientation - Z. P sp ki, M. Storath, D. Sage, M. Unser, "Transforms and Operators for Directional Bioimage Analysis: A Survey," Advances in Anatomy, Embryology and Cell Biology, vol. 219, Focus on Bio-Image Informatics, Springer International Publishing, May 21, 2016.)
*	Requires Canny Edge Detector plugin (https://imagej.nih.gov/ij/plugins/canny/index.html - Tom Gibara)
*	Requires MorphoLibJ v1.6.2 (https://imagej.net/plugins/morpholibj - Legland, D., Arganda-Carreras, I., & Andrey, P. (2016). MorphoLibJ: integrated library and plugins for mathematical morphology with ImageJ. Bioinformatics, 32(22), 3532–3534.
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
 *  Version 2.3.3
 *  - Minor bug fixes
 * 
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Script parameters //
///////////////////////
/*									
#@ String(value = "<html><b> SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS </b> <br>doi:10.1371/journal.pone.0229034</br> <br>https://github.com/oseynnes/SMA</br></html>", visibility="MESSAGE") title

#@ String(value = "NB: The analysis only works with the lower aponeurosis orientated horizontally or downwards to the left", visibility="MESSAGE") text0
#@ Boolean (label="Flip image horizontally", value=false, persist=true, description="Required if image displays proximal side to the right") flip

#@ Boolean (label="Panoramic scan", value=false, persist=true, description="panoramic or regular scan") pano
#@ String (label = "Fascicle model", choices= {"Straight", "Curved_spline", "Curved_circle"}, style="radioButtonHorizontal", description="Assume straight or curved fascicles. Circle method according to Muramatsu et al JAP 2002") geometry
#@ String (label = "Type of analysis", choices= {"Current file", "Open file", "Open folder"}, style="radioButtonHorizontal", description="Analyse single image or several images") analysis
#@ String(label = "Image cropping", choices= {"Automatic (requires identical scan depth)", "Manual"}, style="radioButtonHorizontal", persist=true) cropping
#@ String (choices={"None", "Automatic (requires metadata)", "Manual"}, style="radioButtonHorizontal", persist=true, description="Pixel scaling") scaling

#@ String(value = "<html><br> --------------- Aponeuroses --------------- </br></html>", visibility="MESSAGE") text4	
#@ Integer (label="Tubeness sigma", value=8, style="slider", min=2, max=14, stepSize=2, persist=true, description="Standard deviation of the Gaussian filter. Proportional to aponeurosis thickness") Tsigma
#@ Boolean(label= "Enhance aponeuroses filter", value=false, persist=true, description="run 'Enhance Local Contrast' plugin (CLAHE). Not recommended unless aponeuroses lack contrast") clahe_ap
#@ Integer(label="Length (% of FoV width)", value=80, min=50, max=95, stepSize=5, persist=true, description="Expected length of detected aponeuroses relative to FoV width. default: 80%") apLength
#@ String(label = "Extrapolate from (% of aponeurosis length)", choices = {"100%", "50%"}, persist=true, description="") extrapolate_from
#@ Boolean(label= "Run an aponeurosis detection test (ONLY SINGLE IMAGES)", value=false, persist=true, description="Pause script to try other aponeurosis detection parameters") apo_test

#@ String(value = "<html><br> ---------------- Fascicles ---------------- </br></html>", visibility="MESSAGE") text5
#@ Integer (label="ROI height (% of thickness or 1/2 thickness)", value = 50, style="slider", min=40, max=90, stepSize=5, persist=true) ROIheight
#@ Integer (label="ROI width (% of ROI width)", value = 60, style="slider", min=40, max=90, stepSize=10, persist=true) ROIwidth
#@ Boolean(label= "Enhance fascicle filter", value=false, persist=true, description="run 'Enhance Local Contrast' plugin (CLAHE). Not recommended unless fascicles lack contrast") clahe_fasc
#@ String(label = "Laplacian of Gaussian (sigma)", choices = {"0", "1", "2", "3", "4", "5", "6", "7"}, value=1, persist=true, description="Proportional to fascicle thickness: after running the test, choosing the value yielding the greatest angle is recommended") Osigma

#@ String(value = "<html><br> -------------------- Other -------------------- </br></html>", visibility="MESSAGE") text7
#@ Boolean(label="Print analysis parameters", value=true, persist=true, description="Add analysis parameters to the results table") param
#@ Boolean(label="Display outline of detected fascicle fragments", value=true, persist=true, description="Currently only works with single images") disp_fasc
#@ Boolean(label="Developper mode", value=false, persist=true, description="Disable batch mode and show intermediate steps of analysis") dev
*/

// GLOBAL VARIABLES ------------------------------------------------------------
var inputFile = "";
var output = "";
var input = "";
var extension = "";
var pano_crop = false;
var is_stack = false;
var manual_stack_crop = "Whole stack";
var fixed_ROI = newArray();
var image_titles = newArray();
var slice_n = newArray();
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
var Error = newArray();
var Upper = "";
var Lower = "";
var x_upp = "";
var x_low = "";
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

macro "SMA - Simple Muscle Architecture" {
// description
	plugins = newArray("OrientationJ Measure", "MorphoLibJ...", "Non-local Means Denoising", "Canny Edge Detector");
	check_dependencies(plugins);
	secondary_gui();
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
		
		start = getTime(); // Use this to time how long the whole code takes for one image
		title = getTitle();
		IDraw = getImageID();
		if (is_stack== true) {
			slice = getSliceNumber();
		} else {
			slice = NaN;
		}

		if (flip == true)
			run("Flip Horizontally", "slice");

		if (pano == true) { 
			// ONLY tested with Philips HD11 and Hologic Mach30
			
			if (pano_crop == true) {
				setTool("rotrect");
				if (is("Batch Mode")) {
					setBatchMode(false);
				}
				waitForUser("Select region to crop and press 'OK'");
				run("Duplicate...", " ");
				IDcropped = getImageID();
				selectImage(IDraw); close();
				IDraw = IDcropped;
				if (dev == false) {
					setBatchMode(true);
				}
			}
		}
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Image Cropping - Filter and detect aponeuroses //
		////////////////////////////////////////////////////
		
		counter = 0;
		while (counter == 0) {		
		
			selectImage(IDraw);
//			run("Remove Overlay");
			W = getWidth; 
			H = getHeight;
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");		
		
			
			// cropping
			
			crop_boundaries = perform_cropping(cropping, H, W);
			Le = crop_boundaries[0]; Up = crop_boundaries[1];
	
	
			// aponeuroses detection
	
			ap_filter = filter_aponeuroses(apLength, clahe_ap, Tsigma);
			WFoV = ap_filter[0]; HFoV = ap_filter[1]; IDfilteredAp = ap_filter[2];
	
			detect_aponeuroses(IDfilteredAp, HFoV, Up, Le, analysis, title);
	
			ap_indices = select_indices(extrapolate_from, x_upp, x_low, Upper, Lower);
			upper_idx1 = ap_indices[0]; upper_idx2 = ap_indices[1]; lower_idx1 = ap_indices[2]; lower_idx2 = ap_indices[3];
			
			
			//Upper aponeurosis line (y_1 = b_1 * x + a_1)
			Fit.doFit("Straight Line", Array.slice(x_upp, upper_idx1, upper_idx2), Array.slice(Upper, upper_idx1, upper_idx2));
			a_1 = Fit.p(0); b_1 = Fit.p(1);
			betaUp = atan(b_1)* (180/PI); //angle upper aponeurosis
			//Lower aponeurosis line (y_2 = b_2 * x + a_2)
			Fit.doFit("Straight Line", Array.slice(x_low, lower_idx1, lower_idx2), Array.slice(Lower, lower_idx1, lower_idx2));
			a_2 = Fit.p(0); b_2 = Fit.p(1);
			betaDeep = atan(b_2)* (180/PI); //angle deep aponeurosis
				
			// Plot the curves on top of the existing image	
			selectImage(IDraw);
			applyOverlay(Upper, x_upp);
			applyOverlay(Lower, x_low);
	
			counter = 1;
			// User check of aponeurosis overlay
			if (analysis == "Current file" || analysis == "Open file") {  // Only implemented for single image analysis
				if (apo_test == true && is_stack != true) {  // prevent test on stacks
					setBatchMode(false);
					Dialog.createNonBlocking("Pause");
					Dialog.addChoice("Aponeurosis detection", newArray("Accept and continue", "Try new parameters"), "Accept and continue");
					Dialog.addSlider("Tubeness sigma", 2, 14, Tsigma);
					Dialog.addCheckbox("Enhance aponeuroses filter", clahe_ap);
					Dialog.addSlider("Length (% of FoV width)", 50, 80, apLength);
					Dialog.addChoice("Extrapolate from (% of aponeurosis length)", newArray("50%", "100%"), extrapolate_from);
					Dialog.show();
					status = Dialog.getChoice();
					if (status == "try new parameters") {
						Tsigma = Dialog.getNumber();
						clahe_ap = Dialog.getCheckbox();
						apLength = Dialog.getNumber();
						extrapolate_from = Dialog.getString();
						counter = 0;
						run("Remove Overlay");						
						setBatchMode(true);
					} else {
						setBatchMode(true);  // continue analysis
					}
				}
			}
		}
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
        run("Rotate... ", "  angle=betaDeepR grid=1 interpolation=Bicubic");		

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
			makeRectangle(UAxR[0], maxUAyR+5, (UAxR.length-1)*ROIwidth/100, mindistR*ROIheight*0.01);  // ROI closest to superficial insertion (for curved fascicle mode)
			roiManager("add");
		}
		if (geometry == "Straight") {
			deep_roi_x = UAxR[(UAxR.length-1)*0.5] - (UAxR.length-1)*ROIwidth/100*0.5;  // single ROI is centered along FoV width (for straigth fascicle mode)
		} else {
			deep_roi_x = UAxR[(UAxR.length-1)*(1-ROIwidth/100)];  // ROI closest to deep insertion (for curved fascicle mode)
		}
		makeRectangle(deep_roi_x, minLAyR-mindistR*(ROIheight+5)*0.01, (UAxR.length-1)*ROIwidth/100, mindistR*ROIheight*0.01); //deep			
		roiManager("add");

		// process ROI(s)
		roi_init_n = roiManager("count"); 	

		if (roi_init_n < 1) {
			exit("ROIs for fascicle detection were not defined");
		}
		roi_store = newArray(roi_init_n);										
		for (i=2; i<roi_init_n; i++) {
			selectImage(IDrawR);
			start_n = roiManager("count");
			roiManager("select", i);
			run("Duplicate...", " ");
			IDROI  = getImageID();

			// Preprocessing - enhance structures
			run("32-bit");

			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			run("Subtract Background...", "rolling=50 sliding");
			if (clahe_fasc == true) {
				run("Enhance Local Contrast (CLAHE)", "blocksize=54 histogram=256 maximum=32 mask=*None* fast_(less_accurate)");
			}				
			run("Median...", "radius=2");

			if (Osigma == 0) {
				filt_sigma = 1;  // unspecific value to run the script when the fasccile filter "Osigma" can't be matched
			} else {
				filt_sigma = Osigma;
			}
			run("Tubeness", "sigma=filt_sigma use");
			IDvessel2 = getImageID();
			selectImage(IDROI); close ();
		
			// First directional filter with MorphoLibJ	
			run("Directional Filtering", "type=Max operation=Opening line=20 direction=2");
			horizontal_mask_ID = getImageID();
			if (isOpen("Log")){
				selectWindow("Log"); run("Close");
			}
			imageCalculator("Subtract create", IDvessel2, horizontal_mask_ID);
			filtered_IDROI = getImageID();
			selectImage(IDvessel2); close ();
			selectImage(horizontal_mask_ID); close ();

			// Additional directional filter with FFT
			run("FFT");
			IDFFT2 = getImageID();
			selectImage(IDFFT2);							
			setThreshold(find_lower_thresh(99.5), 255);  // used to be 99.82
			run("Convert to Mask");
			run("Inverse FFT");
			IDinvFFT2 = getImageID();			
			selectImage(IDFFT2); close ();
			selectImage(IDinvFFT2);
			run("8-bit");				
			n_frangi = 1;  // number of times to run the Frangi filter
			for (f = 0; f < n_frangi; f++) {
				temp = getImageID();
				selectImage(temp);
				run("Frangi Vesselness (imglib, experimental)", "number=1 minimum=1.000000 maximum=1.000000");
				selectImage(temp); close();
			}
			IDvessel3  = getImageID();
			selectImage(IDvessel3);	

			// run OrientationJ plugin
			run("Select All");
			run("OrientationJ Measure", "sigma="+Osigma);

			selectImage(IDvessel3); close ();

			// retrieve angle values from OrientationJ log window
			table = getInfo("log"); 					
			lines = split(table, "\n"); 
			headings = split(lines[0], "\t"); 
			values = split(lines[1], "\t");
			selectWindow("Log");
			run("Close");
			roi_store[i] = values[6];
			roi_store[i] = toString(roi_store[i]);
			roi_store[i] = replace(roi_store[i], ",", ".");
			roi_store[i] = parseFloat(roi_store[i]);

			// outline and list fragments of fascicles processed in the analysis
			if (disp_fasc == true && is_stack != true) {  // overlay fragments not available yet during stack analyses
//			if (disp_fasc == true) {  // overlay fragments not available yet during stack analyses
				selectImage(filtered_IDROI);
				run("8-bit");
				run("Auto Local Threshold", "method=Otsu radius=15 parameter_1=0 parameter_2=0 white");
				run("Analyze Particles...", "circularity=0.0-0.7 add");
				// calculate translation offsets
				if (geometry == "Curved_spline" || geometry == "Curved_circle") {
					if (i == 2) {
						IDROI_x = UAxR[0]; IDROI_y = maxUAyR+5;  // superficial ROI
					} else if (i == 3) {
						IDROI_x = deep_roi_x; IDROI_y = minLAyR-mindistR*(ROIheight+5)*0.01;  // deep ROI
					}	
				} else if (geometry == "Straight") {
					IDROI_x = deep_roi_x; IDROI_y = minLAyR-mindistR*(ROIheight+5)*0.01;
				}
				// apply rotation and translation transformations to match original image
				new_n = roiManager("count");
				for (j = start_n; j < new_n; j++) {
				    roiManager("select", j); 
				    getSelectionBounds(x_temp, y_temp, _, _);
				    Roi.move(x_temp + IDROI_x, y_temp + IDROI_y);
				    roiManager("update");
				    roiManager("select", j);
					run("Rotate...", "  angle=betaDeep");
//					run("Rotate...", "rotate angle=betaDeep"); //TODO: fix
					Roi.setStrokeColor("55ffff00");  // set outlines colour to transparent yellow
					roiManager("update");
				}
			}
			selectImage(filtered_IDROI); close ();
		}
			
		selectImage(IDrawR); close ();

		// add outline of fascicle fragments to the final display
		selectImage(IDraw);
		roiManager("select", Array.getSequence(roi_init_n));
		roiManager("delete");
		if (roiManager("count") > 0) {
			run("From ROI Manager");
		}

		roiManager("reset");			

		//////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Calculate fascicle length, pennation angle and thickness //
		//////////////////////////////////////////////////////////////
			
		if (geometry == "Curved_spline" || geometry == "Curved_circle") {
		
			alphaDeep = -roi_store[3]; // with angle due to rotation
			alphaUp = -roi_store[2];		
		    
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


			// centre, radius, curvature and fascicle length
			xs = newArray(Oux,IxFasc,Olx); ys = newArray(Ouy,IyFasc,Oly);
			centre = getCircleCenter(xs, ys);
			x_centre = centre[0]; y_centre = centre[1];
			r = sqrt((x_centre-Oux)*(x_centre-Oux) + (y_centre-Ouy)*(y_centre-Ouy)); //radius
			if (y_centre > 0) { //fascicle curved towards upper aponeurosis
				Crv = 1/r; //curvature
			} else {
				Crv = -1/r;
			}
			
			if (geometry == "Curved_circle") { 
				theta1 = atan2(Ouy - y_centre, Oux - x_centre);
				theta2 = atan2(Oly - y_centre, Olx - x_centre);
				c_x = newArray();
				c_y = newArray();
				for (t = theta1; t <= theta2; t += 0.001) {
				    x = x_centre + r * cos(t);
				    y = y_centre + r * sin(t);
				    c_x = Array.concat(c_x, x);
				    c_y = Array.concat(c_y, y);
				}
				xs = c_x; ys = c_y;  // redefine xs and ys as coordinates of arch						
			} 

			makeSelection("polyline", xs, ys);
			run("Fit Spline");
			run("Interpolate");
			List.setMeasurements;
			Lf = List.getValue("Length");				
			roiManager("add");
			RoiManager.setPosition(slice);
		
		} else { //if fascicles are analysed as straight lines
			alphaDeep = -roi_store[2];			
		    
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
			RoiManager.setPosition(slice);
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
				Wnew=W+left*2;
				run("Canvas Size...", "width=Wnew height=H position=Center zero"); 
				for (i=0; i<n; i++)  {
					roiManager("select", i);
				    roiManager("translate", left, 0);
					roiManager("Update");
				}  
			} else if (Nx2<0 && Nx1>W && left<right) {
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
			slice_n = Array.concat(slice_n, slice);
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
			slice_n = Array.concat(slice_n, slice);
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
			saveAs("tiff", output + file);
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
			slice_n,
			Lower_aponeurosis_orientation,
			Pennation_angle,
			Fascicle_length,
			Curvature,
			Thickness,
			Analysis_duration);
		} else {
			Array.show("Results",
			image_titles,
			slice_n,
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
// FUNCTIONS //
////////////////////

function check_dependencies(deps) {
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
			+"\n- Canny Edge Detector"
			+"\n- MorphoLiJ"
			+"\n- Non-local Means Denoising,"
			+"\n "
			+"\nPlease check that the update sites 'BIG-EPFL', IJPB-plugins, 'Biomedgroup' and 'BioVoxxel' are selected."
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
// Open any file, including DICOM and movies
	splitted = split(File.getName(path), ".");
	if (path.endsWith(".dcm") || splitted.length == 1) {  // extension is .dcm or none
		run("Bio-Formats", "open=path autoscale color_mode=Default open_files rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		unstack_color_channels();
	} else if (path.endsWith(".avi") | path.endsWith(".mp4")) {
		run("Movie (FFMPEG)...", "choose=[path] first_frame=0 last_frame=-1");
	} else {
		open(path);
	}
}


function unstack_color_channels() { 
// convert image stack to RGB image
	title = getTitle();
	getDimensions(im_width, im_height, im_channels, im_slices, im_frames);
	if (im_channels > 1) {
		run("Stack to RGB", "slices");
		run("8-bit");
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
	if (folder_path.length > 0) {  // process folder
		list = getFileList(folder_path);
		open_file(folder_path+File.separator+list[0]);
		if (flip == true) {
			run("Flip Horizontally");
		}
	} else if (is_stack == true) {
		if (fixed_ROI.length == 0) {
			run("Duplicate...", " ");
		} else {  // manual_stack_crop == "Whole stack"
			return fixed_ROI;
		}
	}
	setTool("rectangle"); 
	waitForUser("Select area. Click OK when done");
	Roi.getBounds(x, y, width, height);
	roi_coords = newArray(x, y, width, height);
	if (manual_stack_crop == "Whole stack") { 
		fixed_ROI = Array.concat(fixed_ROI,roi_coords);
	}
//	setBatchMode(true);
	if (folder_path.length > 0 || is_stack == true) {
		close ();
	}
	return roi_coords;		
}


function perform_cropping(cropping, H, W) {
	/**
	 * This function performs cropping on an image based on the input parameters.
	 *
	 * cropping - The type of cropping to be performed; "Manual" or "Automatic".
	 * W - The width of the cropping window; used as an offset in the calculation of the left edge of the cropping rectangle.
	 * H - The height of the cropping window; used as an offset in the calculation of the top edge of the cropping rectangle.
	 *
	 * If the 'cropping' parameter is "Manual", the function performs manual cropping.
	 * If the 'cropping' parameter is "Automatic", the function performs automatic cropping.
	 *
	 * returns an array containing two elements: the left (Le) and top (Up) coordinates of the cropping rectangle.
	 */
    if (cropping == "Manual") {
    	if (dev == false) {
			setBatchMode(false);
		}
        manual_roi = manual_cropping("");
    	if (dev == false) {
			setBatchMode(true);
		}
        x = manual_roi[0]; y = manual_roi[1]; width = manual_roi[2]; height = manual_roi[3];
        Le = x+W*0.005;
        Up = y+H*0.01; // shifts cropped window down by offset
        makeRectangle(Le, Up, width*0.99, height*0.98);
    } else {  // Automatic detection of FoV
        run("Duplicate...", " ");
        IDcopy1  = getImageID();
        run("8-bit");
        run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
        run("Convolve...", "text1=[-1 -1 -1 -1 -1\n-1 -1 -1 -1 -1\n-1 -1 24 -1 -1\n-1 -1 -1 -1 -1\n-1 -1 -1 -1 -1\n] normalize");
        run("Median...", "radius=2");
        run("Auto Local Threshold", "method=Median radius=15 parameter_1=0 parameter_2=0 white");
        run("Options...", "iterations=3 count=1 black do=Close");
        run("Analyze Particles...", "size=10000-Infinity add");
        roiManager("Select", 0);
        getSelectionBounds(x, y, width, height);
        roiManager("delete");
        selectImage(IDcopy1); close();
        run("Select None");
        selectImage(IDraw);
        Le = x+width*0.005;
        Up = y+height*0.01; // shifts cropped window down by offset proportional to height
        makeRectangle(Le, Up, width*0.98, height*0.97);
    }
    return newArray(Le, Up);
}


function enhance_filter(clahe_ap, Tsigma) { 
// Apply various filter to enhance aponeuroses and return image ID of filtered image
	run("8-bit");
    run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
    run("Subtract Background...", "rolling=50 sliding");	
    run("Non-local Means Denoising", "sigma=15 smoothing_factor=1 auto");
    run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 saturate");
    if (clahe_ap == true) {
        run("Enhance Local Contrast (CLAHE)", "blocksize=36 histogram=256 maximum=4 mask=*None* fast_(less_accurate)");  // necessary for faint aponeuroses
    }
    run("Tubeness", "sigma=Tsigma use"); // start with 8
    IDvessel1  = getImageID();
    selectImage(IDvessel1);
    run("8-bit");
    
    return IDvessel1
}


function filter_aponeuroses(apLength, clahe_ap, Tsigma) {
    run("Duplicate...", " ");					
    IDFoV = getImageID();		
    WFoV = getWidth; //Dimensions of the field of view
    HFoV = getHeight;
    setColor(0);
    drawLine(0, 0, WFoV, 0);		
    minLength = apLength/100*WFoV;

    // Preprocessing - enhance structures
	IDvessel1 = enhance_filter(clahe_ap, Tsigma);
	selectImage(IDFoV); close();
	
    // Directional filter
	fragment_length = WFoV / 10;
//	fragment_length = minLength;
	run("Directional Filtering", "type=Max operation=Opening line=fragment_length direction=2");
	selectImage(IDvessel1); close();
	IDvessel1  = getImageID();
	if (isOpen("Log")){
		selectWindow("Log"); run("Close");
	}
    // Preprocessing - Threshold from FFT
    run("FFT");
    IDFFT1 = getImageID();
    selectImage(IDvessel1); close();
    selectImage(IDFFT1);
    Wfft = getWidth;
    Hfft = getHeight;
    setThreshold(find_lower_thresh(99.5), 255);	 // changed to suit MorpholibJ filter
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
    run("Canny Edge Detector", "gaussian=2 low=2.5 high=7.5");	 // changed to suit MorpholibJ filter
    run("Analyze Particles...", "size=0-minLength show=Masks");
    IDmask = getImageID(); // Mask of shorter lines		
    imageCalculator("Subtract create", IDinvFFT1, IDmask);	
    IDfilteredAp = getImageID();
    selectImage(IDinvFFT1); close();
    selectImage(IDmask); close();

    return newArray(WFoV, HFoV, IDfilteredAp);
}



function detect_aponeuroses(IDfilteredAp, HFoV, Up, Le, analysis, title) {
    // Start detecting position of aponeuroses
    selectImage(IDfilteredAp);
    U = round(0.3*HFoV);
    V = U + 1;
    Upper = lineFinder(0, U, 1);  // Call to lineFinder function. The 1 is just a dummy variable so the function can distinguish Upper/Lower
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
    
    if (ok == 1) {
        L = round(0.05*WFoV);
        R = round(0.95*WFoV);
        W2 = R - L;		
        // Adjust aponeuroses to same length
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
    }

}


function select_indices(extrapolate_from, x_upp, x_low, Upper, Lower) {
    // Select indices of aponeuroses lines to extrapolate from
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
//      upper_idx2 = Math.floor(x_upp.length / 2);
//      lower_idx1 = Math.floor(x_low.length / 2);
    }
    return newArray(upper_idx1, upper_idx2, lower_idx1, lower_idx2);
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
	
	// Calculate constants based on the width of field of view (WFoV)
	L = round(0.05*WFoV); // Left boundary
	R = round(0.95*WFoV); // Right boundary
	M = 0.5*WFoV; // Middle point

	// Search for lines from the middle of the image
	makeLine(M, in1, M, in2);
	profile = getProfile();
	
	// Store the indices of potential lines in the 'hold' array
	hold = newArray(1);
	for (j=0; j<profile.length; j++) // First count how many 'lines' there are
		if (profile[j] > 50)  // 50 is a threshold- pixel intensities above 50 are white
			hold = append(hold, j + in1, 0); // Store the indices of each possible line in hold
	hold = Array.slice(hold, 1, hold.length);
	
	// Initialize variables for line tracking
	line = newArray(1); ind1 = 0; ind2 = 1;
	
	// Determine the starting index based on the value of in3 (1 for upper line, 0 for lower)
	if (in3 == 1)
		a = hold.length - 1;
	else
		a = 0;
	
	found = 0;  // Flag to indicate if a suitable line is found
	
	// Iterate until a suitable line is found or the search is exhausted
	while (found != 1 && a >= 0 && a < hold.length) {
		line[0] = hold[a]; // Start with right half of line
		
		// Check just above/below each 'line' for a parallel line
		if(in3 == 1)
			makeLine(M+1, line[0]-21, M+1, line[0]-1);
		else
			makeLine(M+1, line[0]+1, M+1, line[0]+20);
		
		pr = getProfile();
		sum = 0;
		
		// Calculate the sum of pixel intensities in the profile
		for(c=0; c<pr.length; c++)
			sum += pr[c];
		
		if(sum < 50)
			sum += 1; // Do nothing (this is not the right line, try the next one)
		else {
			b = line[0];
			
			// Find the right half of the line
			for(s=M+1; s<R; s++) {
				makeLine(s, b-1, s, b+1); // Assume aponeurosis location varies by no more than +/-1 per pixel
				prof = getProfile();
				
				for(t=0; t<prof.length; t++) {					
					if(prof[t] > 50) {
						b = (b - 1) + t;
						line = append(line, b, 0);  
						t = prof.length; 
						ind2 = s;
					}
				}
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
						ind1 = s;
					}
				} 
			}
		}	
		if (line.length > 2) { // If a suitable line is found, stop the loop
			a = hold.length;
			found = 1;
		}
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
						ind2 = s;
					}
				}
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
						ind1 = s;
					}
				} 
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


function getCircleCenter(xs, ys) {
    // Calculate midpoints of the two lines
    midPoint1x = (xs[0] + xs[1]) / 2;
    midPoint1y = (ys[0] + ys[1]) / 2;
    midPoint2x = (xs[1] + xs[2]) / 2;
    midPoint2y = (ys[1] + ys[2]) / 2;

    // Calculate the slopes of the two lines
    slope1 = (ys[1] - ys[0]) / (xs[1] - xs[0]);
    slope2 = (ys[2] - ys[1]) / (xs[2] - xs[1]);

    // Calculate the slopes of the perpendicular bisectors
    perpBisectorSlope1 = -1 / slope1;
    perpBisectorSlope2 = -1 / slope2;

    // Calculate the y-intercepts of the perpendicular bisectors
    yIntercept1 = midPoint1y - perpBisectorSlope1 * midPoint1x;
    yIntercept2 = midPoint2y - perpBisectorSlope2 * midPoint2x;

    // The center of the circle is the intersection of the two perpendicular bisectors
    x_centre = (yIntercept2 - yIntercept1) / (perpBisectorSlope1 - perpBisectorSlope2);
    y_centre = perpBisectorSlope1 * x_centre + yIntercept1;

    return newArray(x_centre, y_centre);
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
	if (metadata.length > 200) {  // DICOM metadata are typically larger
		if (getInfo("Physical Delta X #1") > 0) {
			scaling_tag = "Physical Delta X #1";
		} else if (getInfo("Physical Delta X") > 0) {
			scaling_tag = "Physical Delta X";
		}
		dcmScale = parseFloat(getInfo(scaling_tag));  // original scale in cm/pxl
		scaling_factor = 1/dcmScale/10;  // scaling factor in pxl/mm
	} else {
		Dialog.create("No metadata");
		Dialog.addRadioButtonGroup("", newArray("No scaling", "Set scale manually", "Exit"), 1, 3, "No scaling")
		Dialog.show();
		action = Dialog.getRadioButton();
		if (action == "No scaling") {
			scaling_factor = 1;
		} else if (action == "Set scale manually") {
			scaling_factor = userscale();  // scaling factor in pxl/mm
		} else {
			exit("Analysis process ended");
		}
	}
	return scaling_factor;  // scaling factor in pxl/mm
}


function userscale() { 
// function description
	run("Select None");
	setTool("line");

	load_settings_file();
	Dialog.createNonBlocking("Manual scaling");
	Dialog.setLocation(0, 0);
	Dialog.addSlider("Scaling distance (cm)", 2, 10, List.getValue("scale_distance"));
	Dialog.addMessage("<html><b> Draw a line over the scaling distance </b> <br> (usually, over the scaling graduation over or on the side of the field of view) </br></html>");
	Dialog.show();
	scale_distance = Dialog.getNumber(); List.set("scale_distance", scale_distance);
	update_settings_file();
	List.clear();

	getLine(x1, y1, x2, y2, lineWidth);
	if (x1 == -1) {
		userscale();
	}
	lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
	return lineLength/(parseInt(scale_distance)*10);  // scaling factor in pxl/mm
}


function handle_stacks() { 
// function description
	list = getList("image.titles");
	if (list.length==0) {
		List.clear();
	    exit("No image windows are open");
	}
	if (nSlices == 1)  // not a stack
		return;
	is_stack = true;
	parent = File.directory;
	sep = File.separator;
	if (endsWith(parent, sep+sep) == true) {
		parent = substring(parent, 0, parent.length-1);
	}
	
	Dialog.createNonBlocking("Image stack");
	Dialog.addCheckbox("Flip image horizontally", flip);
	if (cropping == "Manual") {
		choices1 = newArray("Each frame", "Whole stack");
		Dialog.addRadioButtonGroup("Set manual cropping for:", choices1, 2, 1, "Whole stack");
	}
	choices2 = newArray("Current frame only", "Multiple frames");
	Dialog.addRadioButtonGroup("Process:", choices2, 2, 1, "Current frame only");
	Dialog.addSlider("start", 1, nSlices, 1);
	Dialog.addSlider("end", 1, nSlices, nSlices);
	Dialog.addMessage("Alternatively, enter comma-separated frame numbers");
	Dialog.addString("Discrete frames", "", 4);
	Dialog.addCheckbox("Save processed video (AVI)", false);
	Dialog.addCheckbox("Save results (CSV)", false);
	Dialog.show();
	flip = Dialog.getCheckbox();
	if (cropping == "Manual") {
		manual_stack_crop = Dialog.getRadioButton();
	}
	target = Dialog.getRadioButton();
	start = Dialog.getNumber();
	end = Dialog.getNumber();
	slices_string = Dialog.getString();
	save_to_avi = Dialog.getCheckbox();
	save_to_csv = Dialog.getCheckbox();

	if (target == "Current frame only") {
		start = getSliceNumber(); end = start;
	}
	if (slices_string != "") {
		slices = split(slices_string, ",");
		for (i = 0; i < slices.length; i++) {
			slices[i] = parseInt(slices[i]);
		}
		start = 0; end = slices.length - 1;
	}
	scaleFactor = scale(scaling, "");

	if (dev == false) {
		setBatchMode(true);
	}
	run("Remove Overlay");
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
	Overlay.show;
	if (target != "Current frame only" && save_to_avi == true) {
		avi_path = parent + title + "_processed.avi";
		stack_to_avi(start, end, avi_path);
	}
	if (save_to_csv == true) {
		csv_path = parent + title + "_results.csv";
		saveAs("Results", csv_path);
	}
	if (dev == false) {
		setBatchMode(false);
	}
	exit;
}


function create_settings_log(logpath) { 
// create settings log file if it does not exists
// settings are used for the GUI generated with the Dialog class in the secondary_gui function
	List.clear();
	var_names = newArray(  // list variables names and default values
		"flip", 
		"pano", 
		"pano_crop",
		"geometry", 
		"analysis", 
		"extension",
		"inputFile",
		"input",
		"output",
		"extension",
		"cropping",
		"Tsigma",
		"clahe_ap",
		"apLength",
		"extrapolate_from",
		"ROIheight",
		"ROIwidth",
		"clahe_fasc",
		"n_frangi",
		"Osigma",
		"scaling",
		"scale_distance",
		"param",
		"disp_fasc",
		"dev");
	
	vars = newArray(
		false, 
		false,
		false,
		"Straight", 
		"Current file", 
		".tif",
		getDir("downloads"),
		getDir("downloads"),
		getDir("downloads"),
		".tif",
		"Automatic (requires identical scan depth)",
		10,
		false,
		80,
		"50%",
		"50",
		"60",
		false,
		1,  // number of times to run the Frangi filter
		"1",
		"None",
		"5",
		true,
		true,
		false);
		  
	// pair default settings keys - values in a list
	for (i = 0; i < var_names.length; i++) {
		List.set(var_names[i], vars[i]);
	} 
    f = File.open(logpath);  // command that creates file
	print(f, List.getList);  // add list to file
    File.close(f);
}


function load_settings_file() { 
// function description
	logpath = getDir("plugins")+"SMA_prefs.txt";  // set settings log path
	if(!File.exists(logpath)) {
		create_settings_log(logpath);
	}
	logstring = File.openAsString(logpath);
	logstring_lines = split(logstring, "\n");
	for (i = 0; i < logstring_lines.length; i++) {
		if (logstring_lines[i] > 2) {
			split_line = split(logstring_lines[i], "=");
			List.set(split_line[0], split_line[1]);
		}	
	}
}


function update_settings_file() { 
// update settings
	f = File.open(getDir("plugins")+"SMA_prefs.txt");
	print(f, List.getList);
	File.close(f);
}


function secondary_gui() {
// set additional dialog if required
	if (pano == true || analysis == "Open file" || analysis == "Open folder") {
		load_settings_file();
		Dialog.createNonBlocking("Additional settings");
		
		if (pano == true) {
			Dialog.addMessage("Subselection of panoramic scan");
			Dialog.setInsets(0, 0, 20);
			Dialog.addCheckbox("Pre-crop image", List.getValue("pano_crop"));
		}
		if (analysis == "Open file") {
			Dialog.addFile("Select a file", List.getValue("inputFile"));
		} else if (analysis == "Open folder") {
			Dialog.addDirectory("Input directory", List.getValue("input"));
			Dialog.addDirectory("Output directory", List.getValue("output"));
			Dialog.addChoice("File extension", newArray(".bmp", ".BMP", ".jpg", ".tif", ".png", ".dcm", ""), List.get("extension"));
		}
		
		Dialog.show();
			
		if (pano == true) {
			pano_crop = Dialog.getCheckbox(); List.set("pano_crop", pano_crop);
		}
		if (analysis == "Open file") {
			inputFile = Dialog.getString(); List.set("inputFile", inputFile);
		} else if (analysis == "Open folder") {
			input = Dialog.getString(); List.set("input", input);
			output = Dialog.getString(); List.set("output", output);
			extension = Dialog.getChoice(); List.set("extension", extension);
		}
		
		update_settings_file();
		List.clear();
	}
}


function stack_to_avi(start, end, path) { 
// function description
	run("Flatten", "stack");
	run("Make Substack...", "slices="+start+"-"+end);
	subID = getImageID();
	avi_path = path + "processed.avi";
	run("AVI... ", "compression=JPEG frame=25 save="+path);
	selectImage(subID); close();
}

