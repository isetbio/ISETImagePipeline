cm = cMosaicTreeShrewCreate();

you can get some information for the generated cone mosaic object by typing

properties(cm).

You will see a list of the properties of the cone mosaic. You can get the number of cones in the mosaic via many different properties, for example the coneTypes property (which returns a vector containing the type of each cone in the mosaic, as follows:

conesNum = numel(cm.coneTypes).

Note that the size of activations returned by

cm.compute(theOI, 'nTrials', trialsNum),

will be [trialsNum x timePointsNum x conesNum].

In your simulations, where we do not have a temporal dimension (yet), the returned activations will be [trialsNum x 1 x conesNum].