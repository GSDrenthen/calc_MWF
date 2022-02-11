# calc_MWF

Code used to calculate the myelin-water fraction from 2D GRASE imaging.

* The fourier transform of the excitation pulse is used to determine the imperfect slice profile (by Drenthen et al.: https://doi.org/10.1016/j.neuroimage.2019.04.011).
* The Extended Phase Graph (EPG) is used to estimate the echo intensities (by Weigel: https://doi.org/10.1002/jmri.24619).
* The singular value decomposition filter is used to reduce noise in the multi echo data (by Bydder and Du: https://doi.org/10.1016/j.mri.2006.03.006).
* The Orthogonal Matching Pursuit is used to determine the myelin-water fraction (by Drenthen et al.: https://doi.org/10.1002/mrm.27600)

Contact: g.drenthen@maastrichtuniversity.nl
