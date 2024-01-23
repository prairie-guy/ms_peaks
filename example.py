####################################################################################################
#
#   example.py
#
#   C. Bryan Daniels, cdaniels(at)uchicago.edu
#
#   1/22/2024
#
#   Example: Analysis of Mahdi Mass Spec Data from C -> 4mC Experiments
#
#   Raw data files in `/data`
#
#   Results saved in `/figures`
#
#
####################################################################################################


from fmc import *
from pathlib import Path

data = Path("data")
figures = Path("figures")

Control = read_ms_txt(data/"2024_01_11_C_mer_Control_0_I9_1.txt")
C  = read_ms_txt(data/"2024_01_11_Recipe_C_pH_7.3_K2SO3_950_MeNH2HCl_2024_min_20_0_J8_1.txt")
D1 = read_ms_txt(data/"2024_01_11_10.25_K2S2SO5_666_MAAcid_1400_MeABase_802_min_20_0_K9_1.txt")
D2 = read_ms_txt(data/"2024_01_11_7.77_K2S2SO5_666_MAAcid_1600_MeABase_546_min_20_0_L9_1.txt")
D3 = read_ms_txt(data/"2024_01_11_6.83_K2S2SO5_666_MAAcid_1700_MeABase_418_min_20_0_M9_1.txt")
D4 = read_ms_txt(data/"2024_01_11_7.18_K2S2SO5_666_MAAcid_1661_MeABase_468_min_20_0_N9_1.txt")
D5 = read_ms_txt(data/"2024_01_11_7.14_K2S2SO5_666_MAAcid_1650_MeABase_482_min_20_0_O9_1.txt")
D6 = read_ms_txt(data/"2024_01_11_7.29_K2S2SO5_666_MAAcid_1639_MeABase_496_min_20_0_P9_1.txt")

samples_by_ph = [Control, C, D1, D2, D6, D4, D5, D3]

short_labels = ["Control", "Recipe C", "D1", "D2", "D6", "D4","D5","D3"]

med_labels  = ["C-mer Control",
               "C:  Recipe C pH 7.3 K2SO3 950 MeNH2HCl 2024 min 20",
               "D1: pH 10.25, K2S2SO5 666, MeAAcid 1400, MeABase 802, min 20",
               "D2: pH  7.77, K2S2SO5 666, MeAAcid 1600, MeABase 546, min 20",
               "D6: pH  7.29, K2S2SO5 666, MeAAcid 1639, MeABase 496, min 20",
               "D4: pH  7.18, K2S2SO5 666, MeAAcid 1661, MeABase 468, min 20",
               "D5: pH  7.14, K2S2SO5 666, MeAAcid 1650, MeABase 482, min 20",
               "D3: pH  6.83, K2S2SO5 666, MeAAcid 1700, MeABase 418, min 20"]

full_labels = ["C-mer Control",
               "C: 2024_01_11_Recipe_C_pH_7.3_K2SO3_950_MeNH2HCl_2024_min_20",
               "D1: 2024_01_11_10.25_K2S2SO5_666_MAAcid_1400_MeABase_802_min_20",
               "D2: 2024_01_11_7.77_K2S2SO5_666_MAAcid_1600_MeABase_546_min_20",
               "D6: 2024_01_11_7.29_K2S2SO5_666_MAAcid_1639_MeABase_496_min_20",
               "D4: 2024_01_11_7.18_K2S2SO5_666_MAAcid_1661_MeABase_468_min_20",
               "D5: 2024_01_11_7.14_K2S2SO5_666_MAAcid_1650_MeABase_482_min_20",
               "D3: 2024_01_11_6.83_K2S2SO5_666_MAAcid_1700_MeABase_418_min_20"]

# plot_spec draws the outline of peaks; Good for all peaks
plot_spec(C, figures/"plot_spec_C")

# plot_spke_spec darws filled in peaks; Good for showing fewwer peaks
# call_peaks only returns peaks with intensity > height, 1500 based upon scipy.find_peaks
plot_spike_spec(call_peaks(C, 750), figures/"plot_spike_spec_C")

# plot_multiple peaks shows multiple sample spectrum in one figure
plot_multiple_specs(samples_by_ph, short_labels, figures/"plot_multiple_specs")

# plot_efficiency creates a figure of C -> 4mC efficiency by sample. champ is the prior most efficient sample
plot_efficiency(samples_by_ph, short_labels, figures/"plot_efficiency", champ = "Recipe C")

# table_efficiency creates a table of efficiency by sample
table_efficiency(samples_by_ph, short_labels, figures/"table_efficiency")
