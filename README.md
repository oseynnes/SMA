# SMA #
SMA is an ImageJ macro to automate the analysis of muscle architecture from B-mode ultrasound images. It was developped in colaboration with Dr. Neil Cronin, from the University of Jyväskylä.

A description of the processing operations can be found in the original publication:
https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0229034

Installations steps can be found in the "Supporting information" file published with the article. 
####  NB: for all users, the optional step 3 to add SMA to the *Plugins* menu has changed. Here are updated instructions ####

## Supporting information ##

### Installation instructions for SMA (Simple Muscle Architecture Analysis): An ImageJ/Fiji based macro for automated analysis of B-mode ultrasound images. ###

__Step 1__: Download the Fiji software using the following link. Please read the warning on this page about where to install Fiji, as this can affect your ability to get updates:
https://imagej.net/Fiji/Downloads

__Step 2__: After downloading and running Fiji, you need to gather a few dependencies. First, click Help from the dropdown menu in Fiji, then choose Update. In the Updater window, select Manage update sites. In the window that appears, you will notice that some boxes are already ticked. This means that if/when an update is released for the corresponding programme, your machine will automatically install it at the first opportunity.
Scroll down the list and tick the BIG-EPFL and Biomedgroup update sites. Then, click the Add update site button and scroll to the bottom of the list, where a new update site has been added. Modify the details of the new site as follows (to modify a section, double-click it): Name: SMA
URL: http://sites.imagej.net/SMA Host: webdav:SMA
After entering the details, make sure that the box to the left of the SMA text is ticked. Click Close once you have finished editing. Click Apply changes to confirm the installation. Close Fiji and re-open it. You are now ready to use SMA.

__Optional step 3 (recommended)__: It is convenient to add the SMA macro to the list of available plugins so that you can run SMA from the dropdown menu. Otherwise you’ll need to open the SMA_1_6.ijm file from the Fiji macros folder, and run it manually. **[Updated]** *To add SMA to the Plugins dropdown list start an update via Help > Update. In the ImageJ Updater window, select "Advanced mode". Whith "View all files" as a View option, search for "SMA". If the "Status/Action" column says "Not installed" for the row "plugins/SMA_1_7.ijm", select the row, click "Install" and restart.* **[Updated]** The SMA macro should now appear at the bottom of the dropdown Plugins list.

__Optional step 4__: To analyse movie frames (see limitations in the companion article) with SMA, movies must first be imported into Fiji and converted into image sequences (NB: movies can also be converted in advance using some other software). This is possible if the update site for the FFMEG plugins is added (follow the procedure described above).
A movie file can then be imported by clicking File > Import > Movie (FFMPEG), selecting the movie and accepting the default Import options. Once imported, the movie can be converted into an image sequence in several ways.
One way is to click Image > Stacks > Tools > Make Substack... and then enter the range of frames to be exported. The obtained image stack can then be archived by clicking File > Save as >Image Sequence, and selecting the appropriate file format and saving location.
