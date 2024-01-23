def Peaks(mz, intensity):
    """
    Peaks :: [mz] -> [Intensity] -> Peaks
    Create a Peaks from lists of mz and intensity
    """

def read_ms_mzXML(file):
    """
    read_ms_mzXML :: mzXML_file -> Peaks
    Read mzXML file format and return Peaks
    """

def read_ms_txt(file):
    """
    read_ms_txt :: file -> Peaks
    Read raw ms file text format and returns Peaks
    """

def is_peaks(peaks):
    """
    is_peaks :: Peaks -> Bool
    If peaks is of type Peaks returns True; otherwise returns False
    """
def check_peaks(peaks):
    """
    check_peaks -> Peaks -> IO
    Only returns False if peaks is not of type Peaks; otherwise does nothing
    """

def get_mz(peaks):
    """
    get_mz :: Peaks -> mz
    Returns mz array of peaks
    """

def get_intensity(peaks):
    """
    get_intensity :: Peaks -> Intensity
    Returns intensity of peaks
    """

def empty_peaks(peaks):
    """
    empty_peaks :: Peaks -> Bool
    Returns True if peaks is empty and False if not
    """

def len_peaks(peaks):
    """
    len_peaks :: Peaks -> Int
    Returns length of peaks
    """

def peak(peaks, nth=0):
    """
    peak:: [peaks] -> peaks
    Return nth peak; default is first peak
    """

def mz(peaks, nth=0):
    """
    mz :: peaks -> mz
    Return nth mz; default is first mz
    """

def intensity(peaks, nth=0):
    """
    intensity :: peaks -> Intensity
    Return nth intensity; default is first intensity
    """

def range_peaks(peaks, mz_base, delta=2):
    """
    range_peaks :: Peaks -> Real -> Int -> Peaks
    Return peaks within the range:  (mz_base - delta) < mz < (mz_base + delta)
    """

def area_peak(peaks, mz_base, dist = 0.5):
    """
    area_peak :: Peaks -> mz -> Real -> Int
    Return sum of intensity within range_peaks(peaks, mz_base, dist)
    """

def call_peaks(peaks, height = 1500):
    """
    call_peaks :: Peaks -> Int -> Peaks
    Calls peaks with intensity > height
    """

def group_peaks(peaks, height = 1500, dist=5):
    """
    group_peaks :: Peaks -> Int -> [Peaks]
    Returns a group of peaks, each grouped of which is within a distance dist of each other
    """

def apex_peaks(peaks, height = 1500, dist = 5):
    """
    apexes :: Peaks -> Peaks
    Returns the peak with the highest intensity for each group returned by group_peaks
    Note: Number of apex peaks can very with selection of height and dist. The defaults work well.
          To optimize, consider modifying to apex_peaks to loop through combinations of paramaters.
    """
