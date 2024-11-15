# SMA (Simple Muscle Architecture Analysis)

SMA is an ImageJ macro to automate the analysis of muscle architecture from B-mode ultrasound images. It was developed in collaboration with Dr. Neil Cronin, from the University of Jyväskylä.

A description of the processing operations can be found in the original publication:
*Seynnes, O. R., & Cronin, N. J. (2020). Simple Muscle Architecture Analysis (SMA): An ImageJ macro tool to automate measurements in B-mode ultrasound scans. _PloS One_, _15_(2), e0229034. [https://doi.org/10.1371/journal.pone.0229034](https://doi.org/10.1371/journal.pone.0229034)*

Installations steps can be found in the "Supporting information" file published with the article.  
NB: for all users, the optional step 3 to add SMA to the *Plugins* menu has changed. These instructions are updated bellow.

## Supporting information
Basic user instructions can be found in the original article but there is no real manual for now, except for [this troubleshooting guide](https://oseynnes.github.io/SMA/). Features added after the initial release are only briefly explained in the [changelog](https://github.com/oseynnes/SMA/blob/master/changelog.md). Functions of the starting interface also have a short mouse-over explanation.

### Installation instructions

**__Step 1__**: Download the **Fiji** software using the following link. Please read the warning on this page about where to install Fiji, as this can affect your ability to get updates:
https://imagej.net/Fiji/Downloads

**__Step 2__**: After downloading and running Fiji, you need to gather a few **dependencies**. First, click `Help` from the dropdown menu in Fiji, then choose *Update*. In the Updater window, select `Manage update sites`. In the window that appears, you will notice that some boxes are already ticked. This means that if/when an update is released for the corresponding programme, your machine will automatically install it at the first opportunity.
Scroll down the list and tick the **BIG-EPFL** , **IJPB-plugins**, **BioVoxxel** and **Biomedgroup** update sites.  
If in addition you want to install the **latest version** of SMA (not the **original one**, see next section), click the `Add update` site button and scroll to the bottom of the list, where a new update site has been added. Modify the details of the new site as follows (to modify a section, double-click it): 

Name | URL                         | Host
---- | ----------------------------|-----------
SMA  | https://sites.imagej.net/SMA | webdav:SMA

After entering the details, make sure that the box to the left of the SMA text is ticked.  
Click *Close* once you have finished editing. Click `Apply changes` to confirm the installation. Close Fiji and re-open it. You are now ready to use SMA.

**__Optional step 3 (recommended)__**: It is convenient to **add the SMA macro to the list of available plugins** so that you can run SMA from the dropdown menu. Otherwise you’ll need to open the file from the Fiji macros folder, and run it manually.  

**Latest version:**  
To add the latest version of SMA to the Plugins dropdown list, start an update via `Help > Update`. In the ImageJ Updater window, select `Advanced mode`. With "View all files" as a View option, search for "SMA". If the "Status/Action" column says "Not installed" for the row "plugins/SMA_*version_number*.ijm", select the row, click `Install`. Click `Apply changes` to confirm the installation and restart. The SMA macro should now appear at the bottom of the dropdown Plugins list.  
**Original version (1.7.1):**  
Simply download the file `SMA_171.ijm` and save it in the Fiji `Plugins` folder or in a subfolder. It will appear in the corresponding menu (e.g. Plugins › SMA 171 or Plugins › MyScripts › SMA 171) after restarting the program.

**__Optional step 4__**: **To analyse movie frames** (see limitations in the companion article) with SMA, movies must first be imported into Fiji and converted into image sequences (NB: movies can also be converted in advance using some other software). This is possible if the update site for the **FFMEG** plugins is added (follow the procedure described above).  
**Opening a movie in FijI independently from SMA:**  
A movie file can then be imported by clicking `File > Import > Movie (FFMPEG)`, selecting the movie and accepting the default *Import* options. Once imported, the movie can be converted into an image sequence in several ways.  
One way is to click `Image > Stacks > Tools > Make Substack...` and then enter the range of frames to be exported. The obtained image stack can then be archived by clicking `File > Save as >Image Sequence`, and selecting the appropriate file format and saving location.  
**Opening a movie within SMA** (requires version 2.2):  
Select `Open file` in the `Type of analysis` section. Choose the path from the `Single path file section`.  
