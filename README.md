# Integrated Localization Environment (ILE)

ILE is a software package for segmentation, localization, post-processing and visualization of 2D single molecule localization microscopy (SMLM) data. The package includes post-processing routines for combining signals occuring in consecutive frames, filtering and drift correction. As visualization routines it includes the following methods: scatterplot, histogram binning, Gaussian blurring, triangulation, Voronoi tesselation and local density-based visualization. ILE provides a GUI and allows batch processing. For a quick overview take a look at [ILE's workflow](help/ILE_workflow.png).

## Getting Started

The following sections describe how to install a copy of the software. Detailed instructions can be found in the [manual](help/instructions.pdf).

### Requirements

* Matlab R2012a or newer
	* Statistics and Machine Learning Toolbox 
	* Image Processing Toolbox

* [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin) (a copy is included in this distribution)

### Installing

* download the software package from https://github.com/Jan-NM/ILE
* extract **ILE-master.zip**
* copy the generated **ILE-master** directory into your local Matlab working directory (on Windows machines, usually C:\Users\ "user name"\Documents\MATLAB)
* to use ILE, right click on ILE-master in Matlab's current folder panel, go to "*Add to Path*" and click on "*Selected Folders and Subfolders*"
* type "*startSPDM*" in the command window. Detailed instructions on how to use the software can be found in the [manual](help/instructions.pdf).

## License

ILE is licensed under the GNU GPL - see the [LICENSE](LICENSE) file for details. ILE includes [multiWaitbar](https://de.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar-label-varargin), which comes with a separate license.

## Notes

ILE's segmentation and localization algortihm is based on: Grüll at al., "Accelerating Image Analysis for Localization Microscopy with FPGAs", (2011).

This project is forked from: https://gitlab.com/microscopy/ILE
