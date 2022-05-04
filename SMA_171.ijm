///////////////////////////////////////////////
// SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS //
///////////////////////////////////////////////

//	Requires OrientationJ plugin v16/01/2018 (http://bigwww.epfl.ch/demo/orientation - Z. Püspöki, M. Storath, D. Sage, M. Unser, "Transforms and Operators for Directional Bioimage Analysis: A Survey," Advances in Anatomy, Embryology and Cell Biology, vol. 219, Focus on Bio-Image Informatics, Springer International Publishing, May 21, 2016.)
//	Requires Canny Edge Detector plugin (https://imagej.nih.gov/ij/plugins/canny/index.html - Tom Gibara)
//	Requires Non Local Means Denoise plugin v1.4.6 (https://imagej.net/Non_Local_Means_Denoise - Pascal Behnel, Thorsten Wagner)

//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.

/*
 * 	Change log:
 * 	Version 1.7.1
 * 	- added possibility to extrapolate aponeuroses from 100% or 50% of detected length
 * 	
 * 	Version 1.7
 * 	- added error handling
 * 	- check depencies
 * 	- fixed bug that caused analysis to fail when images are flipped and cropping manual
 */


// OPTION PARAMETERS ------------------------------------------------------------
/*
#@ String(value = "----- SMA - SIMPLE MUSCLE ARCHITECTURE ANALYSIS -----", visibility="MESSAGE") title
#@ String(value = "NB: Analysis only works with lower aponeurosis orientated horizontally or downwards to the left", visibility="MESSAGE") text0
#@ Boolean (label="Flip image horizontally", value=false, persist=false, description="Required if image displays proximal side to the right") flip
#@ String (label = "Type of analysis", choices= {"Image", "Folder"}, persist=false, style="radioButtonHorizontal", description="Analyse single image or several images") analysis
#@ Boolean (label="Print analysis parameters", value=true, persist=true, description="Add analysis parameters to the results table") param
#@ String(value = "------------------ Folder analysis ------------------", visibility="MESSAGE") text1
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File extension", choices = {".bmp", ".BMP", ".jpg", ".tif", ".png"}) extension
#@ String(value = "------------------ Image cropping ------------------", visibility="MESSAGE") text2
#@ String (label = " ", choices= {"Automatic", "Manual"}, style="radioButtonHorizontal", persist=false) cropping
#@ String(value = "--------------- Aponeuroses ---------------", visibility="MESSAGE") text3
#@ Integer (label = "Tubeness sigma (default: 10)", value = 10, description="Standard deviation of the Gaussian filter. Proportional to aponeurosis thickness") Tsigma
#@ String (label = "Extrapolate from (% of detected aponeurosis length)", choices = {"100%", "50%"}, persist=true, description="") extrapolate_from
#@ String(value = "---------------- Fascicles ----------------", visibility="MESSAGE") text4
#@ Integer (label = "Number of ROIs", value = 3, min=2, max=8, persist=false, description="number of ROIs used to detect dominant orientation") ROIn
#@ Integer(label="ROI width (% of FoV width)", value = 60, min=50, max=80, stepSize=10, persist=false) ROIwidth
#@ Integer(label="ROI height (% of thickness)", value = 90, min=40, max=90, stepSize=10, persist=false) ROIheight
#@ Boolean (label="Automatic thresholding", value=true, persist=true, description="Based on % of pixels above threshold in the power spectrum") autThresh
#@ Integer (label = "Manual threshold value", value = 175, persist=false) manThresh
#@ String (label = "Laplacian of Gaussian (sigma)", choices = {"0", "1", "2","3", "4", "5", "6", "7", "Test"}, value=4, persist=false, description="Proportional to fascicle thickness: after running the test, choosing the value yielding the greatest angle is recommended") Osigma
#@ String (label = "Main orientation method", choices = {"max", "median", "mean (not recommended)"}, persist=true, description="Choose value retained from analyses of ROIs ('max' is recommanded in most cases)") Pa
#@ String(value = "---------------- Pixel scaling ----------------", visibility="MESSAGE") text5
#@ Boolean (label="Scaled measurements (requires identical scan depth)", value=false, persist=false, description="Requires scaling bar and identical scan depth when analysing multiple images") scaling
#@ Integer(label="Scan depth (cm)", min=3, max=7, style="scroll bar", value = 5) depth
*/

// GLOBAL VARIABLES ------------------------------------------------------------
var Scan = newArray();
var ALPHA = newArray();
var Lower_aponeurosis_orientation = newArray();
var Pennation_angle = newArray();
var Fascicle_length = newArray();
var Thickness = newArray();
var Analysis_duration = newArray();
var Image_cropping = newArray();
var Tubeness_sigma = newArray();
var ROI_n = newArray();
var ROI_w = newArray();
var ROI_h = newArray();
var Thresholding = newArray();
var OrientationJ_sigma = newArray();
var Orientation = newArray();
var Error = newArray();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This macro processes an opened image or all the images in a folder and any subfolders.
macro "SMA - Simple Muscle Architecture Analysis" {

	missing = newArray();
	deps = newArray("OrientationJ Measure", "Non-local Means Denoising", "Canny Edge Detector");
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
	
	if (isOpen("Results")) { //Close results window (obtained from previous analysis, should be kept out of loop for stacks)
		selectWindow("Results");
	    run("Close");
	}
	if (isOpen("ROI Manager")) {
		selectWindow("ROI Manager");
	 	run("Close");
		}

	if (analysis == "Image") {
		title = getTitle();
		run("Select None");
		run("Remove Overlay");
		if (cropping == "Manual") {
			setTool("rectangle");
			waitForUser("Select area. Click OK when done");
			Roi.getBounds(x, y, width, height);
			if (scaling == true) {
				run("Select None");
				setTool("line");
				waitForUser("Select scaling line. Click OK when done");
				getLine(x1, y1, x2, y2, lineWidth);
				lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
				scaleFactor = lineLength/(depth*10);
			}
			setBatchMode(true);
			singleImageAnalysis();
		} else {
			if (scaling == true) {
				run("Select None");
				setTool("line");
				waitForUser("Select scaling line. Click OK when done");
				getLine(x1, y1, x2, y2, lineWidth);
				lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
				scaleFactor = lineLength/(depth*10);
			}
			setBatchMode(true);
			singleImageAnalysis();
		}
	} else {
		count = 0;
		countFiles(input);
		n = 0;
		n2 = 1;
		if (count == 0)
			exit("This folder does not contain any "+extension+" file.");
		if (cropping == "Manual") {
			list = getFileList(input);
			open(input+File.separator+list[0]);
			if (flip == true)
				run("Flip Horizontally");
			setTool("rectangle");
			waitForUser("Select area. Click OK when done");
			Roi.getBounds(x, y, width, height);
			if (scaling == true) {
				run("Select None");
				setTool("line");
				waitForUser("Select scaling line. Click OK when done");
				getLine(x1, y1, x2, y2, lineWidth);
				lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
				scaleFactor = lineLength/(depth*10);
			}
			close ();
			setBatchMode(true);
			processFolder(input);
		} else {
			if (scaling == true) {
				run("Select None");
				setTool("line");
				waitForUser("Select scaling line. Click OK when done");
				getLine(x1, y1, x2, y2, lineWidth);
				lineLength = sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
				scaleFactor = lineLength/(depth*10);
			}
			setBatchMode(true);
			processFolder(input);
		}
		selectWindow("Processed");
		run("Close");
	}


	function countFiles(input) {
	  list = getFileList(input);
	  for (i=0; i<list.length; i++) {
		if (endsWith(list[i], extension))
	    	count++;
	  }
	}

	function processFolder(input) {
		list = getFileList(input);
		list = Array.sort(list);
		run("Text Window...", "name=Processed");
		for (i = 0; i < list.length; i++) {
			if (endsWith(list[i], extension)) {
				print("[Processed]", list[i] +" (scan #"+n2++ +" out of " +count+ ")" +"\n");
				if (File.isDirectory(input + File.separator + list[i]))
					processFolder(input + File.separator + list[i]);
				if (endsWith(list[i], extension))
					processFile(input, output, list[i]);
			}
		}
	}

	function processFile(input, output, file) {
		open(input+File.separator+file);
		singleImageAnalysis();
	}

	//---------------------------------------------------------
	function singleImageAnalysis() {
		//Start of analysis
		if (flip == true)
			run("Flip Horizontally");
		start = getTime(); // Use this to time how long the whole code takes for one image
		title = getTitle();
		IDraw = getImageID();
		W = getWidth;
		H = getHeight;
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		
		//---------------------------------------------------------
		// Image Cropping
		if (cropping == "Manual") {
			Le = x+W*0.005;
			Up = y+H*0.01; // shifts cropped window down by 1% of its height
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
			selectImage(IDcopy1);
			close ();
			run( "Select None" );
			selectImage(IDraw);
			close("\\Others");
			Le = x+width*0.005;
			Up = y+height*0.01; // shifts cropped window down by 1% of its height
			makeRectangle(Le, Up, width*0.98, height*0.97);
		}

		//---------------------------------------------------------
		//Filter and detect aponeuroses
		run("Duplicate...", " ");
		IDFoV = getImageID();
		WFoV = getWidth; //Dimensions of the field of view
		HFoV = getHeight;
		minLength = 0.8*WFoV;
		run("8-bit");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		run("Subtract Background...", "rolling=50 sliding");
		run("Non-local Means Denoising", "sigma=15 smoothing_factor=1 auto");
		run("Bandpass Filter...", "filter_large=40 filter_small=3 suppress=None tolerance=5 saturate");
		run("Enhance Local Contrast (CLAHE)", "blocksize=36 histogram=256 maximum=4 mask=*None* fast_(less_accurate)");
		run("Tubeness", "sigma=Tsigma use"); // currently recommended: 10
		IDvessel1  = getImageID();
		selectImage(IDvessel1);
		run("8-bit");
		run("FFT");
		IDFFT1 = getImageID();
		selectImage(IDvessel1);
		close ();
		selectImage(IDFFT1);
		Wfft = getWidth;
		Hfft = getHeight;

		// automatic threshold
		percentage = 99.9; // percentage of FFT pixels to threshold out
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
		setThreshold(lowThresh, 255);
		setOption("BlackBackground", true);
		run("Convert to Mask");
		setColor(0);
		makeRectangle(Wfft/2+2, 0, Wfft/2, Hfft/2);
		fill();
		makeRectangle(0, Hfft/2+1, Wfft/2-1, Hfft/2);
		fill();
		run("Inverse FFT");
		IDinvFFT1 = getImageID();
		selectImage(IDFFT1);
		close ();
		selectImage(IDinvFFT1);
		run("Canny Edge Detector", "gaussian=2 low=1 high=7.5");
		run("Analyze Particles...", "size=0-minLength show=Masks");
		IDmask = getImageID(); // Mask of shorter lines
		imageCalculator("Subtract create", IDinvFFT1,IDmask);

		// Start detecting position of aponeuroses
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

		err = "";
		if (Upper.length < 10){
			ok = 0;
			if (analysis == "Image") {
				exit("Could not detect upper aponeurosis. Please try again with different value of tubeness sigma or try manual cropping");
			} else {
				err = "Could not detect upper aponeurosis. Please try again with different value of tubeness sigma or try manual cropping";
				Error = Array.concat(Error, title +":"+err);
				continue;
			}
		} else if (Lower.length < 10) {
			ok = 0;
			if (analysis == "Image") {
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
				Lower = Array.copy(Lower2); }
			else if (Upper.length < Lower.length) {
				maxL = Lower.length;
				Upper2 = newArray(Lower.length);
				for (t=0; t<Upper.length; t++)
					Upper2[t] = Upper[t];
				i1 = Upper[0];
				i2 = Upper[Upper.length-1];
				islope = (i2 - i1) / Upper.length;
				for (p=Upper.length; p<Lower.length; p++)
					Upper2[p] = Upper2[p-1] + islope;
				Upper = Array.copy(Upper2); }
			else
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

			selectImage(IDraw);
			applyOverlay(Upper, x_upp); // Plot the curves on top of the existing image
			applyOverlay(Lower, x_low);
			close("\\Others");

			
			if (extrapolate_from == "100%") {
				upper_idx2 = x_upp.length-1;
				lower_idx1 = 0;
			} else {
				upper_idx2 = x_upp.length/2-1;
				lower_idx1 = x_low.length/2-1;
			}
			//Upper aponeurosis line (y_1 = b_1 * x + a_1)
			b_1 = (Upper[0] - Upper[upper_idx2]) / ( x_upp[0] - x_upp[upper_idx2] ); //slope based on 1/2 aponeurosis
			a_1 = Upper[0] - b_1 * x_upp[0];
			//Lower aponeurosis line (y_2 = b_2 * x + a_2)
			b_2 = (Lower[lower_idx1] - Lower[x_low.length-1]) / ( x_low[lower_idx1] - x_low[x_low.length-1] ); //slope based on 1/2 aponeurosis
			a_2 = Lower[0] - b_2 * x_low[0];
			beta = atan(b_2)* (180/PI); //angle

			//Rename aponeuroses coordinates arrays
			UAx = Array.copy(x_upp); //upper aponeurosis x values
			UAy = Array.copy(Upper); //upper aponeurosis y values
			LAx = Array.copy(x_low); //lower aponeurosis x values
			LAy = Array.copy(Lower); //lower aponeurosis y values

			//---------------------------------------------------------
			// Measure fascicle orientation

			Array.getStatistics(UAx, minUAx, maxUAx, meanUAx, stdUAx); //statistics of upper aponeurosis array
			Array.getStatistics(UAy, minUAy, maxUAy, meanUAy, stdUAy);
			Array.getStatistics(LAx, minLAx, maxLAx, meanLAx, stdLAx); //statistics of lower aponeurosis array
			Array.getStatistics(LAy, minLAy, maxLAy, meanLAy, stdLAy);
			mindist = minLAy-maxUAy; //smallest distance between aponeuroses
			maxdist = maxLAy-minUAy; //largest distance between aponeuroses

			// Set parameters needed to draw ROIs to detect dominant orientation
			par1 = newArray(ROIn); // x coordinate
			par2 = newArray(ROIn); // y coordinate
			height = newArray(ROIn); // Height of ROI
			boxW = ROIwidth/100 * maxL; // ROI width is proportional to aponeuroses overlay
			ref = newArray(ROIn);
			for (h=0; h<ROIn; h++) {
				ref[h] = (maxL-boxW)*h/(ROIn-1);
				cut = 1-ROIheight/100;
				mindTop = (LAy[L + ref[h]] - UAy[L + ref[h]]) * cut; // 0.1 is the proportion of top thickness left out
				if (b_2 <= 0) { //Assumption: the upmost right point on the lower aponeurosis has coordinates (maxLAx, minLAy)
					mindBot =  LAy[ref[h]] - LAy[boxW-1 + ref[h]];
				} else {
					mindBot =  LAy[boxW-1 + ref[h]] - LAy[ref[h]];
				}
				par1[h] = L + ref[h] + Le;
				par2[h] = UAy[L + ref[h]] + mindTop;
				height[h] = (LAy[L + ref[h]] - UAy[L + ref[h]]) - 1.1*mindTop - mindBot; //height of ROI, problem in some cases (eg lower ap. orientated downwards to the right)

				makeRectangle(par1[h], par2[h], boxW, height[h]);
				roiManager("add");
			}

			n = 0;
			cont = 0;
			n = roiManager("count");
			if (n > 0) {
				roiManager("Show All");
				cont = 1; }
			if (cont > 0) {
				store = newArray(ROIn);
				for (i=0; i<n; i++)  {
					roiManager("select", i);
					run("Duplicate...", " ");
					IDROI  = getImageID();
					run("32-bit");
					run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
					run("Subtract Background...", "rolling=50 sliding");
					run("Non-local Means Denoising", "sigma=15 smoothing_factor=1 auto");
					run("Median...", "radius=2");
					getStatistics(area, mean, min, max, std, histogram);
					run("FFT");
					IDFFT2 = getImageID();
					selectImage(IDFFT2);
					if (autThresh == true) {
						percentage = 99.82; // percentage of FFT pixels to threshold out
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
					} else {
						lowThresh = manThresh;
					}
					setThreshold(lowThresh, 255);
					run("Convert to Mask");
					run("Inverse FFT");
					IDinvFFT2 = getImageID();
					selectImage(IDFFT2);
					close ();
					selectImage(IDinvFFT2);
					run("8-bit");
					run("Tubeness", "sigma=1 use");
					IDvessel2  = getImageID();
					selectImage(IDvessel2);
					run("Select All");
					if (Osigma == "0")
						run("OrientationJ Measure", "sigma=0.0");
					else if (Osigma == "1")
						run("OrientationJ Measure", "sigma=1.0");
					else if (Osigma == "2")
						run("OrientationJ Measure", "sigma=2.0");
					else if (Osigma == "3")
						run("OrientationJ Measure", "sigma=3.0");
					else if (Osigma == "4")
						run("OrientationJ Measure", "sigma=4.0");
					else if (Osigma == "5")
						run("OrientationJ Measure", "sigma=5.0");
					else if (Osigma == "6")
						run("OrientationJ Measure", "sigma=6.0");
					else if (Osigma == "7")
						run("OrientationJ Measure", "sigma=7.0");
					else {
						run("OrientationJ Measure", "sigma=0.0");
						run("OrientationJ Measure", "sigma=1.0");
						run("OrientationJ Measure", "sigma=2.0");
						run("OrientationJ Measure", "sigma=3.0");
						run("OrientationJ Measure", "sigma=4.0");
						run("OrientationJ Measure", "sigma=5.0");
						run("OrientationJ Measure", "sigma=6.0");
						run("OrientationJ Measure", "sigma=7.0");
						exit();
					}

					selectImage(IDvessel2);
					close ();
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
					selectImage(IDinvFFT2);
					close ();
					selectImage(IDROI);
					close ();
				}
			}

			Array.getStatistics(store, min, max, mean, stdDev);
			store1 = Array.sort(store);
			if (store1.length % 2 == 1) //if n ROI is odd integer
				median = store1[(store1.length-1)/2];
			else
				median = (store1[(store1.length+1)/2] + store1[store1.length/2-1]) / 2;
			if (Pa == "median")
				alpha = median;
			else if (Pa == "max")
				alpha = min;
			else
				alpha = mean;


			//---------------------------------------------------------
			// Calculate fascicle length, pennation angle and thickness

			//Pennation angle
			theta = (-alpha) - beta;

			//Composite fascicle line  (y_3 = b_3 * x + a_3)
			b_3 = tan(-alpha * (PI / 180));
			if (b_2 <= 0) {
				a_3 = minLAy - b_3 * maxLAx; //Assumption: the upmost right point on the lower aponeurosis has coordinates (maxLAx, minLAy)
			} else {
				a_3 = maxLAy - b_3 * maxLAx;
			}

			//Fascicle-upper aponeurosis intersection coordinates
			Ix = (a_1 - a_3) / (b_3 - b_1);
			//yy = b_3 * xx + a_3
			Iy = a_1 + b_1 * Ix;

			//Offset coordinates on lower aponeurosis
			Hx1 = (Ix + b_2 * Iy - b_2 * a_2) / ((b_2 * b_2) + 1);
			Hy1 = b_2 * ((Ix + b_2 * Iy - b_2 * a_2) / ((b_2 * b_2) + 1)) + a_2;
			Olx = maxLAx + ((minLAx - Hx1)/2);
			if (b_2 <= 0) {
				Oly = minLAy + ((maxLAy - Hy1)/2);
			} else {
				Oly = maxLAy + ((minLAy - Hy1)/2);
			}


			//Offset composite fascicle line  (y_3 = b_3 * x + a_3)
			b_o = tan(-alpha * (PI / 180));
			a_o = Oly - b_o * Olx;

			//Offset coordinates on upper aponeurosis
			Oux = (a_1 - a_o) / (b_o - b_1);
			Ouy = a_1 + b_1 * Oux;

			//unscaled fascicle length
			Lf = sqrt((Olx-Oux)*(Olx-Oux) + (Oly-Ouy)*(Oly-Ouy));
			makeLine(Olx, Oly, Oux, Ouy);
			roiManager("Add");
			n = roiManager("count");

			//adjust fascicle and ROIs display
			Omx = (Olx+Oux)/2; //x coordinate of mid fascicle
			mFoVx = UAx[(UAx.length)/2]; //x coordinate of mid FoV
			if (Oux<0 || Olx>W) {
				roiManager("select", n-1);
				roiManager("translate", -(Omx-mFoVx), 0);
				roiManager("Update");
				getLine(Nx1, Ny1, Nx2, Ny2, NlineWidth); // right most start point is Nx2,Ny2
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

			//Fascicle length
			if (scaling == true)
				Lf = Lf / scaleFactor;

			//Thickness (temporary).
			sta = minOf(x_upp.length, x_low.length);
			Th = newArray(sta);
			sum = 0;
			for (i=0; i<Th.length; i++){
				Th[i] = Lower[i] - Upper[i];
				sum += Th[i];
			}

			if (scaling == true) {
				Th = (sum/Th.length)/scaleFactor;
			} else {
				Th = sum/Th.length;
			}
		}

		// Analysis duration
		Duration = (getTime() - start) / 1000;

		//---------------------------------------------------------
		// Save variables and figures

		// Save variables of interest
		Scan = Array.concat(Scan, title);
		ALPHA = Array.concat(ALPHA, alpha);
		Lower_aponeurosis_orientation = Array.concat(Lower_aponeurosis_orientation, beta);
		Pennation_angle = Array.concat(Pennation_angle, theta);
		Fascicle_length = Array.concat(Fascicle_length, Lf);
		Thickness = Array.concat(Thickness, Th);
		Analysis_duration = Array.concat(Analysis_duration, Duration);
		Image_cropping = Array.concat(Image_cropping, cropping);
		Tubeness_sigma = Array.concat(Tubeness_sigma, Tsigma);
		ROI_n = Array.concat(ROI_n, ROIn);
		ROI_w = Array.concat(ROI_w, ROIwidth);
		ROI_h = Array.concat(ROI_h, ROIheight);
		if (autThresh == true) {
			Thresholding = Array.concat(Thresholding, "auto");
		} else {
			Thresholding = Array.concat(Thresholding, manThresh);
		}
		OrientationJ_sigma = Array.concat(OrientationJ_sigma, Osigma);
		Orientation = Array.concat(Orientation, Pa);

        // Save figures with overlays
        if (analysis == "Folder") {
			selectImage(IDraw);
			run("From ROI Manager");
			run("Flatten", "stack");
			saveAs(extension, output + File.separator + file);
			roiManager("reset");
			close("*");
		}
		// Convert 0 values to NaN for final display
		ALPHA = arrTidy(ALPHA);
		Lower_aponeurosis_orientation = arrTidy(Lower_aponeurosis_orientation);
		Pennation_angle = arrTidy(Pennation_angle);
		Fascicle_length = arrTidy(Fascicle_length);
		Thickness = arrTidy(Thickness);
	}

	//---------------------------------------------------------
	// Display results
	if (param == false) {
		Array.show("Results",
		Scan,
		Lower_aponeurosis_orientation,
		Pennation_angle,
		Fascicle_length,
		Thickness,
		Analysis_duration);
	} else {
		Array.show("Results",
		Scan,
		Lower_aponeurosis_orientation,
		Pennation_angle,
		Fascicle_length,
		Thickness,
		Analysis_duration,
		Image_cropping,
		Tubeness_sigma,
		ROI_n,
		ROI_w,
		ROI_h,
		Thresholding,
		OrientationJ_sigma,
		Orientation);
	}

    if (analysis == "Folder") {
    	list = getFileList(output);
		displayStack(output);
    }
	function displayStack(output) {
		stack1=output+File.separator;
		run("Image Sequence...", "open=&stack1 sort use");
	}

	if (Error.length > 0) {
		Array.show("Images not processed", Error);
	}

	setBatchMode(false);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS ---------------------------------------------------------------

function arrTidy(Arr) {
	for (m=0; m<Arr.length; m++) {
		if (Arr[m] == 0)
			Arr[m] = NaN;
	}
	return Arr;
}

//////////

//in1: vert start point of line, in2: vert end point of line, in3: 1 for upper line, 0 for lower
function lineFinder(in1, in2, in3) {
	//U = round(0.3*HFoV); V = U + 1;
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

//////////

function compDiff() {
oz = newArray(p.length-1);
for(t=0; t<p.length-1; t++){
	oz[t] = p[t+1] - p[t]; // Compute diff pixel by pixel
	if(oz[t] != 0) oz[t] = 1; // Convert to binary
} return oz; }

//////////

function applyOverlay(cInput, x) {
	run("Select All");
	makeSelection("polyline", x, cInput);
	run("Fit Spline");
	run("Overlay Options...", "stroke=green width=4");
	run("Add Selection...");
	run("Select None"); }

//////////
// arr = array to append to, value = value to append, place = 1 for start, 0 for end
function append(arr, value, place) {
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
