####################################################################################################
#
#   fmc.py
#
#   Analysis of Mass Spectrum Data in Study of C -> 4mC (in DNA oligo 5'-A-G-C-G-A-3')
#
#
#   ms_peaks.py provides general primitives to work with data type: Peaks
#
#   This library is specific to for functions associated with C -> 4mC
#
#   C. Bryan Daniels, cdaniels@nandor.net
#
#   1/22/2024
#
####################################################################################################


####################################################################################################
#
# Description of primary data type: Peaks
#
# Peaks are a tuple of two numpy arrays, of which the first array contains values of mz and the
# second the values of intensity. The data type is not strictly enforced as a class,
# but could easily be implemented as such. (Types are loosely defined in Haskell Style)
#
# Type Peaks = ([mz], [intensity])
#
# Type mz = [Float]
#
# Type mz = [Int]
#
####################################################################################################


from ms_peaks import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#def max_peak(peaks, mz_base, delta=2):

def c_peak(peaks,delta=2):
    """
    c_peak :: Peaks -> Peaks
    Return C apex Peaks (mw = 1508)
    """
    return max_peak(peaks, 1508, delta)

def fmc_peak(peaks,delta=2):
    """
    fmc_peak :: Peaks -> Peaks
    Return 4mC apex Peaks (mw = 1522)
    """
    return max_peak(peaks,1522, delta)

def fec_peak(peaks,delta=2):
    """
    fmc_peak :: Peaks -> Peaks
    Return 4ethylC apex Peaks (mw = 1536)
    """
    return max_peak(peaks,1536, delta)

def efficiency(fn_start, fn_end, peaks):
    """
    efficiency :: Function -> Function -> Peaks
    """
    try:
        i_start = intensity(fn_start(peaks))
        i_end   = intensity(fn_end(peaks))
        return abs(round(i_end/(i_start + i_end),2))
    except:
        return 0

def efficiency_fmc(peaks): return efficiency(c_peak, fmc_peak, peaks)

def efficiency_fec(peaks): return efficiency(c_peak, fec_peak, peaks)

def plot_multiple_specs(peaks_list, labels_list, title=None, save_as = None, xrange=(1500, 1600), figsize=(20, 12), annotate=True):
    n = len(peaks_list)
    fig, axs = plt.subplots(n, 1, figsize=figsize, squeeze=False)
    #if title: plt.title(title) # For multiple figures, title is not working
    for i, peaks in enumerate(peaks_list):
        mz, intensity = peaks
        ax = axs[i, 0]
        ax.vlines(mz, ymin=0, ymax=intensity)
        ax.set_xlim(*xrange)
        ax.grid(color='gray', linestyle='-', linewidth=0.5)
        #ax.set_xlabel("m/z")
        #ax.set_ylabel("Intensity")
        title = labels_list[i]
        ax.text(0.5, 0.8, title, horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12, color='black')
        if annotate:
            for (x,y) in zip(get_mz(apex_peaks(peaks)), get_intensity(apex_peaks(peaks))):
                ax.annotate(round(x,1), xy=(x,y/2), textcoords="offset points", xytext=(-5,10), ha='right')
    plt.tight_layout()
    if save_as:
        plt.savefig(f"{save_as}.png")
        plt.close(fig)

def plot_spec(peaks, title = None, save_as=None, xrange = (1500,1600),  figsize = (20,10), annotate = True):
    if len_peaks(peaks) < 500:
        return("`plot_deep_spec` is for original ms data. Use `plot_spec` instead")
    mz, intensity = peaks
    plt.figure(figsize = figsize)
    plt.plot(mz, intensity)
    plt.xlim(*xrange)
    plt.grid(color='gray', linestyle='-', linewidth=0.5)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    if title: plt.title(title)
    if annotate:
        for (x,y) in zip(get_mz(apex_peaks(peaks)), get_intensity(apex_peaks(peaks))):
            plt.annotate(x, xy=(x,y/2), textcoords="offset points", xytext=(-5,10), ha='right')
    if save_as:
        plt.savefig(f"{save_as}.png")
        plt.close()

def plot_spike_spec(peaks, title=None, save_as=None,xrange = (1500,1600),  figsize = (20,10), annotate = True):
    mz, intensity = peaks
    plt.figure(figsize = figsize)
    plt.vlines(mz, ymin=0, ymax=intensity)
    plt.xlim(*xrange)
    plt.grid(color='gray', linestyle='-', linewidth=0.5)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    if title: plt.title(title)
    if annotate:
        for (x,y) in zip(get_mz(apex_peaks(peaks)), get_intensity(apex_peaks(peaks))):
            plt.annotate(x, xy=(x,y), textcoords="offset points", xytext=(-5,10), ha='right')
    if save_as:
        plt.savefig(f"{save_as}.png")
        plt.close()

def table_efficiency(peaks_list, labels_list, fname, figsize=(5,4)):
    data = list(zip(labels_list, map(efficiency, peaks_list)))
    columns = ('Sample', 'Efficiency')
    rows = range(1, len(data) + 1)
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('tight')
    ax.axis('off')
    ax.table(cellText=data, colLabels=columns, loc='center')
    plt.savefig(f"{fname}.png")
    plt.close(fig)

    
def plot_efficiency(peaks_list, labels_list, fname, champ = None, figsize = (6,4)):
    if len(peaks_list) != len(labels_list):
        return "error: The length of 'peaks' and 'labels' must be the same."
    if not champ is None and not champ in labels_list: return(f"error: champ = {champ} not in {labels_list}")
    if not champ is None:
        idx = list(range(len(labels_list))) if champ is None else labels_list.index(champ)
        champ_value = round(efficiency(peaks_list[idx])*100, 0)
        labels_list = labels_list[:idx] + labels_list[idx+1:]
        peaks_list = peaks_list[:idx] + peaks_list[idx+1:]

    plt.rcParams.update({'font.size': 10})
    plt.figure(figsize=figsize)
    if not champ is None:
        plt.axhline(y= champ_value, color='r', linestyle='-', label = champ)
        plt.legend()
    for peaks, label in zip(peaks_list, labels_list):
        eff = round(efficiency(peaks)*100,0)
        plt.vlines(x = label, ymin=0, ymax = eff, linewidth = 5.0 )
    plt.xlabel("Samples")
    plt.ylabel("Efficiency (%)")
    plt.title("Efficiency of C -> 4mC Conversion")
    plt.ylim(0, 100)
    plt.savefig(fname)
    plt.close()




####################################################################################################
#
# Testing
#
####################################################################################################


def label_peaks(peaks):
    ap = apex_peaks(peaks)
    c = c_peak(peaks) # Call c_peak(peaks) NOT c_peaks(apex); group_peaks is too restrictive
    fmc = fmc_peak(peaks) # same

    mz_ap, mz_c, mz_fmc = np.round(get_mz(ap)), round(mz(c)), round(mz(fmc))

    idx_c = [mz in harmonics(mz_c, "K") for mz in mz_ap]
    idx_fmc = [mz in harmonics(mz_fmc, "K") for mz in mz_ap]

    mz_ap = get_mz(ap)
    print(f"Apex Peaks: {mz_ap}" )
    print(f"C Peak: {mz(c)}")
    print(f"C Potassium Peaks: {mz_ap[idx_c]}")

    print(f"4mC Peak: {mz(fmc)}")
    print(f"4mC Potassium Peaks: {mz_ap[idx_fmc]}")

def harmonics(mz, element="K", harmonic = None):
    hs = {"K":38, "Na":22} # Note: harmonic = mw - 1
    if harmonic is None: harmonic = hs[element]
    return round(mz,0) + harmonic * np.arange(5)[1:]

# def efficiency(peaks, height = 1500, dist = 5):
#     c, fc = c_peak(peaks, height, dist), fmc_peak(peaks, height, dist),

#     if len_peaks(c) != 1: return(f"{c} is not peak 1508")
#     c = intensity(c)

#     if len_peaks(fc) > 1: return(f"{c} is not peak 1522")
#     fc = 0 if empty_peaks(fc) else intensity(fc)



def scan_ubisulfide(peaks, delta = 5, dist = 1):
    c = c_peak(peaks)
    fmc = fmc_peak(peaks)
    for d in range(83 - delta, 83 +delta +1) :
        print(d, round(np.sum(get_intensity(range_peaks(C, mz(c) + d, dist)))/np.sum(get_intensity(range_peaks(C, mz(fmc), dist))),3))

def test_labels():
    for label, sample in zip(ph_labels_short[1:], ph_order[1:]):
        print(label)
        label_peaks(sample)
        print("++++++++++++++++++++++++++++++++++++")
