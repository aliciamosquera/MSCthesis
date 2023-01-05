# Analysis of dynamic three-dimensional images for sperm motion studies - An investigation of boar spermatozoa’s kinematics
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

### Images available
Three different datasets were provided by the laboratory in Barcelona. One was recorded in December 2019 (dataset0 equivalent to 9_2) and the other two in December 2021 (dataset1 equivalent to 00068 and dataset2 equivalent to 00064). Dataset 0 was mainly used to get acquainted with the data and see what possibilities we had with it, hence the name. It presented several problems, such as fluorescence in certain background sections, a flow in the whole sample and an imaging artefact similar to a mirroring effect. Since these problems were fixed in datasets 1 and 2, those were the ones used
to analyse the motion. The parameters of the three datasets are specified in [datasets.json](data/datasets.json)

### Processing done outside of this code

- To estimate the noise in the datasets we calculated the intensity’s standard
deviation in several homogeneous regions (both background and foreground) with ImageJ, resulting in the documents [noise_std_...](data)
- The datasets were too big to compute altogether (2.5 GB each), so we divided
them into smaller snippets with ImageJ. Since our aim was to study the
motion of individual cells, it made sense to fragment the datasets into snippets
of cells. The idea was that each section depicted the whole time series of one
cell. In practice, our samples were densely packed with sperm cells, so it was
not possible to define a volume that contained the whole motion of only one
cell. Nonetheless, we prioritised having the entire time series of one cell in the
same movie and later segmenting the different cells. 14 snippets made the cut from datasets 1 & 2.
- Segmentation of the head was done manually off of one snippet to have a baseline. Given that
we could not perform a 3D segmentation ourselves, we segmented three slices for each coordinate plane (xy, yz, zx). Each slice was chosen at random and accepted if we could see a cell in it with the naked eye. For this we made use of [Photopea](https://www.photopea.com/).
- Random forest classifier with ImageJ to segment all the snippets.

## Automatic cell tracking system
The pipeline was created based on one fragment of dataset 1, although it was
later generalised to be able to compute the same parameters for any other
fragment of either dataset 1 or 2.

The outline of the image processing goes as follows:
1. Image segmentation to obtain the contour and centroid position of each sperm head
    We applied three different techniques:

    1.1. one consisting of Gaussian scale-space smoothing and later thresholding

    1.2. random forest classifier with the tool Labkit from ImageJ

    1.3. manual segmentation to compare the former two

2. Linking the position of the same cell in subsequent image frames to obtain the trajectory of each cell
3. Calculating sperm motility parameters from the obtained trajectories 