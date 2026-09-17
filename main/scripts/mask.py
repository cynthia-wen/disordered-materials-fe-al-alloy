import shutil
import numpy as np
import h5py

# file paths
INTENSITY = 'fitting_js/structures_20/relaxed_08/average_diffuse/symm.h5'
MULTIPLIER = 'reciprocal_space_multiplier.h5'
OUTPUT = 'fitting_js/structures_20/relaxed_08/average_diffuse/masked_diffuse_intensity.h5'

# reading data
intensity = h5py.File(INTENSITY, 'r')['/data'][:]
multiplier = h5py.File(MULTIPLIER, 'r')['/data'][:]

if intensity.shape != multiplier.shape:
    raise ValueError(
        f'shape mismatch: intensity {intensity.shape}, '
        f'multiplier {multiplier.shape}.'
    )

# mask
keep = multiplier != 0
n_masked = int((~keep).sum())
print(f'masked       : {n_masked} of {keep.size} voxels '
      f'({100 * n_masked / keep.size:.1f}%)')

masked = np.where(keep, intensity.astype(np.float64), np.nan) # intensity is kept where multiplier is non-zero, otherwise NaN

finite = masked[np.isfinite(masked)]
if finite.size:
    print(f'kept range   : {finite.min():.4g} to {finite.max():.4g}')

# save output file (original grid metadata is preserved)
shutil.copyfile(INTENSITY, OUTPUT)

with h5py.File(OUTPUT, 'r+') as f:
    if f['/data'].dtype == masked.dtype:
        f['/data'][...] = masked
    else:
        # integer or float32 source cannot hold float64 NaN 
        del f['/data']
        f.create_dataset('data', data=masked)

print(f'wrote        : {OUTPUT}')