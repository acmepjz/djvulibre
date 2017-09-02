Fork of DjVuLibre
=================

This is a fork of DjVuLibre, minidjvu, etc. at SourceForge.

The original README of DjVuLibre is [here](README), and the original README of minidjvu is [here](minidjvu/README).

What's new compared to original DjVuLibre and minidjvu:
* switched to CMake build system
* added FreeImage support
* added file list options to cjb2, djvm, minidjvu, etc.
* added some options and experimental features to cjb2:
  * encoding using existing dictionary
  * multipage encoding and generating dictionary
* added some options to cpaldjvu: lossy option, specify palette, etc.

Additional utilities:
* tsv2djvu: convert .tsv generated by Tesseract to .dsed.

Compile
-------

Use CMake to generate makefile.

Optional libraries:
* FreeImage
* FLANN (not yet)
* Leptonica (not yet)
* OpenCV (not yet)
* eigen3 (not yet)