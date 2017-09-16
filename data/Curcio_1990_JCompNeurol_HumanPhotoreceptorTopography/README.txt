Data retrieved from:

	https://info.cs.uab.edu/sloan/PRtopo/

Last checked 9/6/2017

Data from this paper:

Curcio, Christine A., et al. "Human photoreceptor topography." Journal of comparative neurology 292.4 (1990): 497-523.

The file 'curcioConeDensityPerSqMm.mat' contains the mean cone density values (and associated support in mm eccentricity) from these data.

Curcio's measurements were made on the retinal surface in units of mm and counts per mm^2. While the source data provides an eccentricity support vector in units of degrees, this differs a bit from the values produced by our mm --> deg conversion function (which in turn is derived from Watson's paper). So, we load the data in mm / mm^2 and then convert within the routine.