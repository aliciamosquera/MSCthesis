# Analysis of dynamic three-dimensional images for sperm motion studies
## An investigation of boar spermatozoa’s kinematics
Alicia Mosquera Rodríguez

May 20, 2022

Supervisors: Prof. Jon Sporring and Asst. Prof. Amin Doostmohammadi

MSc in Physics, University of Copenhagen

## Abstract
Negative trends in fertility and quality of sperm have been reported, with evidence
to suspect that endocrine disrupters impact sperm motion. In particular,
they can provoke alternations between progressive and hyperactivated motility,
lowering the chances of fertilisation. We monitored boar sperm cells to
determine their motion descriptors and correlate them to the motility regime.
Our raw data were dynamic three-dimensional images of boar spermatozoa
obtained by light sheet fluorescence microscopy. We describe the methodology
followed to go from images, with the mere information being pixel intensities,
to physical measures. We defined a centroid and basis for the spermatozoon’s
head. Then, we computed motion descriptors such as the angular rotation
and speed by following the head’s movement through time. We found that
the head does not have an intrinsic rotation while swimming, but it shows
symmetric bending about its three coordinate axes. The imaging technique
and segmentation enlarged the actual head, resulting in average head dimensions
five times the value given in the literature and high uncertainty in
speed determination. Therefore, we could not confidently assure the motility
regime of the cells in our sample, but we saw that they move ballistically even
without a concentration gradient to guide them. In conclusion, this thesis conceived
a methodology to measure the sperm cell’s speed, rotation, and motion
regime from 3D videos. At the same time, the imaging technique fulfilled the
expectations of being an effective tool for studying sperm dynamics.

## Technical details

### Data
Three different datasets were provided by the laboratory in Barcelona. One
was recorded in December 2019 (dataset 0) and the other two in December
2021 (datasets 1 & 2). Dataset 0 was mainly used to get acquainted with the
data and see what possibilities we had with it, hence the name. It presented
several problems, such as fluorescence in certain background sections, a flow
in the whole sample and an imaging artefact similar to a mirroring effect.
Since these problems were fixed in datasets 1 and 2, those were the ones used
to analyse the motion. The parameters of the three datasets are specified in _datasets.json_.
- The folder **data** contains the processing developed with dataset 0
- The folder **newdata** contains the processing developed with datasets 1 and 2

###