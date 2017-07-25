Active Shape Models for face detection
======
![](/Media/Faces_MultiResolution_horizontal.png "Variations in the gray-level model")

This project was part of my work for *Advanced Digital Image Processing* at the University of Iowa during spring 2017. You can read my [final report for the class here](/Media/ADIP_ActiveShapeModels_FinalReport.pdf). The report was a tad rushed, my apologies!

![](/Media/Video/ASM_FaceDetection_24-Jul-2017_MUCT.gif "Finding a face using the MUCT layout")

### Usage ###
After cloning this repository, run the [Example_FindFace](Example_FindFace.m) script for a walkthrough demonstration of how to use this ASM code for locating a face in an example image.

### Background ###
Here is the original [Cootes et al. paper.](http://www.sciencedirect.com/science/article/pii/S1077314285710041) PDFs of the paper are available elsewhere online if you don't have access to the journal. Here is a link to the [faces training set](http://robotics.csie.ncku.edu.tw/Databases/FaceDetect_PoseEstimate.htm#Our_Database_) I annotated to train my model.

Manipulating the weights on the 1st and 2nd principal components deforms the face shape within an allowable range of variation.

<img src="/Media/Faces_PC_Variations.png" width="800">
