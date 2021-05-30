Attracted Growth Cone Search (AGCS)

This script is used to support the addition of biologically realistic jitter to curves based on histological data. The basic principle of this method is to utilized hypothetical growth cones to control the meander angles which are added in the new line. This is preferable to random jitter which can add extreme meander angles and overlapping segments if not properly controlled. AGCS can be better used to control the parametric distribution of meander angles without disrupting the distribution of segment lengths beyond that which is reasonable.

![AGCS Example](https://github.com/bingsome/AttractedGrowthConeSearch/blob/master/docs/animation.gif)

![AGCS comparison with Random Jitter](https://github.com/bingsome/AttractedGrowthConeSearch/blob/master/docs/Meanders.png)

Originally this package was written to aid the construction of biologically realistic neural models for use in computational studies of extracellular electrical stimulation but it can be used for diverse neural modeling problems.

The core script is found in AttractedGrowthConeSearch.py and is executed thusly:
```

AttractedGrowthConeSearch() takes no arguments but should instatiated as a class and then the .jitter_streamline() method should be used with the following arguments to execute the algorithm:

line = list of control points that you wish to jitter with AGCS
angle_kde = a kernel density estimate of angles which AGCS can sample to create its growth cones
mode=this is the mode of operation of the AGCS algorithm which can be either 'forward' or 'backward'. 'backward' is the default argument and is much more stable, though forward is more biologically realistic mechanism.

The method returns a list of new control points to which AGCS has been applied.
```

Example code can be found in /tests.

