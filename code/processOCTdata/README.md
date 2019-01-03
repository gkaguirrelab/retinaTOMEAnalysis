We have collected OCT data using a Heidelberg Spectralis system. A horizontal and vertical macular scan was collected for the left and right eye for each subject. The data were exported from the OCT system and saved in `.E2E` format (which is proprietary to Heidelberg) and in `.vol` format.

The `.vol` files were copied to a separate location and then processed by an operator (Kara Cloud) using OCT Explorer v5.0 (on a Macintosh). This analysis yields (for each `.vol` file) a file named `_Surfaces_Retina-JEI-Final.xml`.

The routine `xmlConversionWrapper.m` is used to convert the `.xml` file to a `.mat` file. The wrapper makes use of the routine `xml2volmask.m` which is present within the `octExplorerSupport` toolbox, and was written by Jin Gahm from the LONI group, USC.

The resulting `.mat` file has the dimensions 768x97x496, corresponding to the vertical, axial (depth) and horizontal diimensions. Each voxel is given an integer value from 0 - 11, corresponding to a retinal layer (with a value of zero indicating that the voxel does not reside within the retina). The depth dimension is in units of mm, while the transverse (horizontal and vertical) dimensions are in units of degrees of visual angle; our macular acquisitions were 30° wide.
