Use active shape models from Cootes et al. to locate faces in images.

Here is the original [Cootes et al. paper.](http://www.sciencedirect.com/science/article/pii/S1077314285710041) PDFs of the paper are available elsewhere online if you don't have access to the journal.

Here is a link to the [faces training set](http://robotics.csie.ncku.edu.tw/Databases/FaceDetect_PoseEstimate.htm#Our_Database_) I annonated to train my model.

Manipulating the weights on the 1st and 2nd principal components deforms the face shape within the allowable range of variation.
![Variation of the 1st and 2nd face principal components](/Media/Faces_PC_Variations.png)