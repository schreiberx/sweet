#!/usr/bin/env python3

import sys
import numpy as np

from mule.postprocessing.JobData import JobData
from mule.postprocessing.SphereDataSpectral import SphereDataSpectral
import mule.postprocessing.SphereDataOperators as SphereDataOperators

debug_active = False




job_dir = sys.argv[1]
phi_file = sys.argv[2]
vrt_file = sys.argv[3]
div_file = sys.argv[4]

# Load job's data
job_data = JobData(job_dir)
job_data_flattened = job_data.get_flattened_data()

# TODO: Determine this automagically
rsphere = 6371220
grav = 9.81

# Load fields
phi_spec = SphereDataSpectral(phi_file, setup_physical=False)
vrt_spec = SphereDataSpectral(vrt_file, setup_physical=False)
div_spec = SphereDataSpectral(div_file, setup_physical=False)

# Setup transformations
sh = SphereDataOperators.SphereDataOperators(rsphere=rsphere)
sh.setup(phi_spec.file_info, anti_aliasing=False)

# Setup transformations
sh_aa = SphereDataOperators.SphereDataOperators(rsphere=rsphere)
sh_aa.setup(phi_spec.file_info, anti_aliasing=True)



######################################################
# Compute
# Ke = 1/2 * m * V^2
######################################################

# Compute u and v, prepared for anti-aliasing
u_phys_data, v_phys_data = sh_aa.vrtdiv2uv(vrt_spec.data_spectral, div_spec.data_spectral)

# Compute
#    u*u + v*v
# and apply anti-aliasing
V2_phys_data = sh.spec2phys(sh_aa.phys2spec(u_phys_data*u_phys_data + v_phys_data*v_phys_data))

# Get mass (which we relate to the height of the SWE)
# m = geopot. / g
m_phys_data = sh.spec2phys(phi_spec.data_spectral)/grav

# Finish computation of
# Ke = 1/2 * m * V^2
# in physical space
ke = 0.5 * m_phys_data * V2_phys_data



######################################################
# Determine spectrum
######################################################



def dist_bucket(real_m, verbose=False):
    """
    Compute the bucket and distribution of mode real_m which is not integer

    Now we need to split things up into buckets, 0th mode bucket

    m = 0 ... 0.5 ... 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ...
        |      |               |               |
        | bck0 |    bucket1    |    bucket2    |   bucket3
        |      |               |               |

    Change the real_m to:

    m = 0.5 .. 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ... 3.5 ...
         |      |               |               |
         | bck0 |    bucket1    |    bucket2    |   bucket3
         |      |               |               |
    """

    if verbose:
        print(" +++ using real mode ", real_m)

    real_mh = real_m + 0.5

    # Now things are easy...

    bucket_a_num = int(real_mh - 0.5)
    bucket_b_num = int(real_mh + 0.5)

    bucket_a_weight = 1.0 - (real_mh - 0.5 - bucket_a_num)
    bucket_b_weight = real_mh + 0.5 - bucket_b_num

    if verbose:
        print(" +++ bucket ", bucket_a_num, " gets ", bucket_a_weight)
        print(" +++ bucket ", bucket_b_num, " gets ", bucket_b_weight)

    assert(bucket_a_num >= 0)
    assert(bucket_b_num == bucket_a_num + 1)
    assert(bucket_a_weight >= 0 and bucket_a_weight <= 1.0)
    assert np.allclose(bucket_a_weight + bucket_b_weight, 1.0)
    return bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight





def dist_bucket_array(real_m, verbose=False):
    """
    Compute the bucket and distribution of mode real_m which is not integer

    Now we need to split things up into buckets, 0th mode bucket

    m = 0 ... 0.5 ... 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ...
        |      |               |               |
        | bck0 |    bucket1    |    bucket2    |   bucket3
        |      |               |               |

    Change the real_m to:

    m = 0.5 .. 1.0 ... 1.5 ... 2.0 ... 2.5 ... 3.0 ... 3.5 ...
         |      |               |               |
         | bck0 |    bucket1    |    bucket2    |   bucket3
         |      |               |               |
    """

    if verbose:
        print(" +++ using real mode ", real_m)

    real_mh = real_m + 0.5

    # Now things are easy...

    bucket_a_num = np.array(real_mh - 0.5, dtype=int)
    bucket_b_num = np.array(real_mh + 0.5, dtype=int)

    bucket_a_weight = 1.0 - (real_mh - 0.5 - bucket_a_num)
    bucket_b_weight = real_mh + 0.5 - bucket_b_num

    if verbose:
        print(" +++ bucket ", bucket_a_num, " gets ", bucket_a_weight)
        print(" +++ bucket ", bucket_b_num, " gets ", bucket_b_weight)

    assert(np.greater_equal(bucket_a_num, 0).all())
    assert(np.equal(bucket_b_num, bucket_a_num+1).all())

    assert(np.greater_equal(bucket_a_weight, 0).all())
    assert(np.less_equal(bucket_a_weight, 1).all())

    assert np.allclose(bucket_a_weight + bucket_b_weight, 1.0)
    return bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight


if debug_active:
    bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket(0, True)
    assert(bucket_a_num == 0)
    assert(bucket_b_num == 1)
    assert np.allclose(bucket_a_weight, 1.0)
    assert np.allclose(bucket_b_weight, 0.0)

    bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket(0.1, True)
    assert(bucket_a_num == 0)
    assert(bucket_b_num == 1)
    assert np.allclose(bucket_a_weight, 0.9)
    assert np.allclose(bucket_b_weight, 0.1)

    bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket(2.1, True)
    assert(bucket_a_num == 2)
    assert(bucket_b_num == 3)
    assert np.allclose(bucket_a_weight, 0.9)
    assert np.allclose(bucket_b_weight, 0.1)

    bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket(2.9, True)
    assert(bucket_a_num == 2)
    assert(bucket_b_num == 3)
    assert np.allclose(bucket_a_weight, 0.1)
    assert np.allclose(bucket_b_weight, 0.9)


    bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket_array(np.array([0, 0.1, 2.1, 2.9]), True)
    assert np.equal(bucket_a_num, [0, 0, 2, 2]).all()
    assert np.equal(bucket_b_num, [1, 1, 3, 3]).all()
    assert np.allclose(bucket_a_weight, [1.0, 0.9, 0.9, 0.1])
    assert np.allclose(bucket_b_weight, [0.0, 0.1, 0.1, 0.9])



print("Resolution in physical space: ", ke.shape)
ke_spec = np.fft.rfft(ke, axis=1)
print("Resolution after longitudinal FT transformation: ", ke_spec.shape)

def spec_to_buckets_iter(mode_numbers, modes_data, buckets):
    """
    Iterate over all Fourier modes
    m here relates to the number of waves
    """
    for m in range(len(mode_numbers)):
        # Compute real mode (including shortening by being closer to poles)
        # There would be additional number of waves, hence we need to divide by this
        bucket_a_num, bucket_a_weight, bucket_b_num, bucket_b_weight = dist_bucket(mode_numbers[m])

        ampl = np.abs(modes_data[m])

        if bucket_a_num < len(buckets):
            buckets[bucket_a_num] += ampl

            if bucket_b_num < len(buckets):
                buckets[bucket_b_num] += ampl



def spec_to_buckets_fast(mode_numbers, modes_data, buckets):

    bucket_a_num_, bucket_a_weight_, bucket_b_num_, bucket_b_weight_ = dist_bucket_array(mode_numbers)
    ampl_ = np.abs(modes_data)

    for m in range(len(buckets)):
        if bucket_a_num_[m] < len(buckets):
            buckets[bucket_a_num_[m]] += ampl_[m]

            if bucket_b_num_[m] < len(buckets):
                buckets[bucket_b_num_[m]] += ampl_[m]



num_buckets = ke_spec.shape[1]-1
print("Setting up ", num_buckets, "spectral buckets")
buckets = np.zeros(num_buckets)

# Iterate over all longitude stripes
for i in range(ke_spec.shape[0]):
    print("Lat: ", sh.lats[i])

    # Compute scalar to multiply modal number with
    # scaling factor \in [
    s = np.cos(sh.lats[i])

    m_ = range(num_buckets)
    real_m_ = np.array(m_) / s

    if debug_active:
        a = np.zeros_like(buckets)
        spec_to_buckets_iter(real_m_, ke_spec[i], a)

        b = np.zeros_like(buckets)
        spec_to_buckets_fast(real_m_, ke_spec[i], b)

        print(a-b)
        assert np.allclose(a, b)


    if 1:
        spec_to_buckets_fast(real_m_, ke_spec[i], buckets)


#
# Plot results
#

import matplotlib.pyplot as plt

def modes_to_wavelengths(modes):
    return np.pi * 2.0 * rsphere / modes


def plot_buckets(buckets, label):
    # Compute wavelengths
    _ = np.arange(len(buckets))
    _[0] = -1.0 # avoid div/0
    wavelengths = modes_to_wavelengths(_)

    plt.plot(wavelengths[1:], buckets[1:], label=label)

    plt.gca().invert_xaxis()
    plt.title("Kinetic Energy spectrum")
    plt.xlabel("Wavelength")
    plt.ylabel("Amplitude")

    plt.xscale("log")
    plt.yscale("log")


def plot_k_lines(ax):

    k_ = np.array([100, 200])
    wl_ = modes_to_wavelengths(k_)

    #
    # k^-3
    #
    s = 1e13
    y_ = np.power(k_, -3.0)*s
    line = plt.plot(wl_, y_, linestyle="solid", color="gray")
    x = line[0].get_xdata()[1]
    y = line[0].get_ydata()[1]

    ax.annotate(
        "k^-3",
        xy=(x * 1.05, y * 0.35),
        color=line[0].get_color(),
        size=10,
    )

    #
    # k^-(5/3)
    #
    s *= np.power(k_[0], -3.0)/np.power(k_[0], -5.0/3.0)
    y_ = np.power(k_, -5.0/3.0)*s
    line = plt.plot(wl_, y_, linestyle="solid", color="gray")
    x = line[0].get_xdata()[1]
    y = line[0].get_ydata()[1]

    ax.annotate(
        "k^-5/3",
        xy=(x * 1.05, y * 1.35),
        color=line[0].get_color(),
        size=10,
    )


fig, ax = plt.subplots(figsize=(6, 4))

plot_buckets(buckets, "test")
plot_k_lines(ax)
plt.legend()

plt.tight_layout()
plt.savefig("output.pdf")
