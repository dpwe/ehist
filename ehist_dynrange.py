#!/usr/bin/env python
"""
ehist_dynrange.py

Calculate the dynamic range of a sound in different frequency bands 
as a way of detecting noisy, non-phone-channel sound examples.

2014-04-09 Dan Ellis dpwe@ee.columbia.edu
"""

import numpy as np
# For filter
import scipy.signal

def frame(x, window, hop):
    """ Convert vector x into an array of successive window-point segments, 
        stepped by hop
        Done with stride_tricks, no copying, but have to lose the final 
        part-window 
    """

    nframes = 1 + int( np.floor( (len(x)-window)/hop ))
    shape = (nframes, window)
    strides = (x.strides[0] * hop, x.strides[0])
    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)

def stftm(signal, nfft, window, hop):
    """ calculate the short-time fourier transform magnitude """
    frames = frame(signal, window, hop)

    # apply frame window to each frame
    window = np.hanning(window)
    wframes = frames * window

    return np.abs(np.fft.rfft(wframes, int(nfft)))

def percentile(data, pcntl):
    """
    v = percentile(d,n)
    Return for each col of v the n'th percentile where 0<n<1
    2004-10-04 dpwe@ee.columbia.edu
    """
    nr = np.size(data, axis=0)

    x = np.sort(data, axis=0)
    return x[int(np.floor(pcntl*nr)),]

def sadmask(dBD, pcntl=0.1, marg=10.0, collar=10):
    """
    M = sadmask(dBD, pcntl, marg, collar)
    Build a sad mask from a dB spectrogram
    pcntl is the point on the CDF to choose the low
    threshold (default 0.1 = 10th percentile).
    margin is the dB margin above the low threshold to count as the
    speech active threshold (dflt 10.0).
    collar is the number of frames to spread the speech active
    regions in both directions, and also median filter window 
    (dflt 10). 
    2013-05-27, 2014-01-03 Dan Ellis dpwe@ee.columbia.edu
    """
    # total energy is taken as mean across all dB histo bins (!)
    mdBD = np.mean(dBD, axis=0)

    thresh = percentile(mdBD, pcntl)
    M1 = mdBD > (thresh + marg)

    # grow it out this far in both directions
    #collar = 10   # with 16 ms hop, this is +/- 160ms
    M = M1
    for i in range(collar):
        # OR together versions of the mask shifted both ways
        M = M | np.r_[False, M[:-1]] | np.r_[M[1:],False]

    # Median filter too
    #M = medianf(M,collar);
    M = (scipy.signal.medfilt(M, collar)>0)
    # I'm pretty sure this does nothing if applied after spreading

    # but in all cases exclude frames with supernegative bins (digital
    # silence)
    M = M & (np.median(dBD, axis=0)>0)
    
    return M

def histc(data, edges):
    """
    counts = histc(data, edges)
    forms histograms of the values in each row of data
    bins lists a single set of (monotonic) band edges
    counts[i] returns the number of values falling into
    the bin defined by edges(i-1) < val <= edges(i)
    edges[-1] is taken to be -Inf, and edges(end) = Inf
    2014-04-09 Dan Ellis dpwe@ee.columbia.edu
    """
    rows = np.size(data, axis=0)
    edges = np.r_[edges, np.inf]
    counts = np.zeros( (rows, len(edges)) )
    lastcount = np.zeros( rows );
    for ix, edge in enumerate(edges):
        newcount = np.sum(data <= edge, axis=1)
        counts[:,ix] = newcount - lastcount
        lastcount = newcount
    return counts

def ehist_calc(d, sr, domask=True, offset=80.0, rng=100.0):
    """
    ehist_calc - calculate energy histogram for a waveform
    [X,DNM,CV,F] = ehist_calc(N, domask, offset, rng, melbins)
    X return a 2D histogram of levels/dB vs. freq for an entire
    conversation side. 
    offset is a dB offset added to every sgram value (dflt 70).
    rng is the largest dB value (scalar), or the range of dB
    values to use (pair), or a set of dB bin centers (vector).
    DNM returns the parent spectrogram, with noise gating if selected
    Originally called uttftr.m .
    CV returns a covariance matrix.
    F returns the frequencies for each bin.
    2013-05-25 Dan Ellis dpwe@ee.columbia.edu
    """

    dlen = len(d)

    dbbinwidth = 1.0
    dbbins = range(int(np.round(rng/dbbinwidth)))

    # Choose fft size; 1024 for sr = 16000
    nfft = int(np.power(2.0, 
                        np.round(np.log(1024.0 * sr/16000.0)/np.log(2.0))))
    nhop = nfft/2
    nsgframes = 1+ int(np.floor((dlen-nfft)/nhop))
    fr = sr/nhop
    # Signal spectrogram
    SG = stftm(d, nfft, nfft, nhop).T
    # lose the nyquist bin as it messes up the histogram
    SG = SG[:nfft/2,]

    DN = 20.0*np.log10(SG)
    # Actual frequency bin values
    F = [float(x) / nfft * sr for x in range(nfft/2)]

    # collapse to 0..80 dB range (from top)
    DN = DN + offset

    if domask:
        # Build a sad mask from a dB spectrogram
        percentl = 0.1
        marg     = 10.0
        collar   = 5
        # not sure why I have to take the [0] row of output of np.nonzero
        DN = DN[:, np.nonzero(sadmask(DN, percentl, marg, collar))[0]]

    # Calculate subband energy histogram in each row, in 1 dB bins
    ncols = np.size(DN, axis=0)
    X = np.zeros( (len(dbbins)+1, ncols), np.float)
    #for i in range(ncols):
    #    X[:,i] = histc(DN[i,], dbbins)
    X = histc(DN, dbbins).T

    # Covariance matrix
    #CV = np.cov(np.exp(DNM/8.68588963806504))

    return X, F

def histpercentile(data, pctls=0.5):
    """
    P = histpercentile(D,T)
    Columns of D are histograms.  Return the 100*T'th percentile
    (default 0.5) for each column in P.
    If T is a vector, return multiple rows, one for each percentile.
    2013-06-15 Dan Ellis dpwe@ee.columbia.edu
    """

    nrows, ncols = np.shape(data)
    npctls = len(pctls)
    P = np.zeros( (npctls, ncols), np.float)

    # Just one cumsum
    cs = np.cumsum(data, axis=0) / np.tile(np.sum(data, axis=0), (nrows,1))
    for tx in range(npctls):
        P[tx,] = np.sum(cs <= pctls[tx], axis=0)

    return P

########### main function ##############

def ehist_dynrange(d, sr):
    """
    D = ehist_dynrange(d, sr)
    Calculate an energy histogram for the specified soundfile
    and report the dB difference between 5th and 95th percentiles
    for each of the 5 octaves 
    125-250 Hz, 250-500 Hz, 500-1 kHz, 1-2 kHz, 2-4 kHz.
    2014-04-09 Dan Ellis dpwe@ee.columbia.edu
    """

    domask = 1
    offset = 90.0
    rng = 120.0
    X, F = ehist_calc(d, sr, domask, offset, rng)

    # Calculate difference of 5th and 95th percentiles
    pctls = [0.05, 0.95]
    hp = histpercentile(X, pctls)
    dynrng = np.diff(hp, axis=0)[0]

    # Find bin edges
    fmin = 125.0
    fmax = 4000.0
    fthis = fmin
    nix = int(np.floor(1.0 + np.log(fmax/fmin)/np.log(2.0)))
    ix = np.zeros( nix, np.int )
    for i in range(nix):
        ix[i] = np.amax(np.nonzero(np.less_equal(F, fthis))[0])
        fthis *= 2.0

    # Calculate average dyn ranges in those fbin ranges
    Dout = np.zeros(nix-1, np.float)
    for i in range(nix-1):
        Dout[i] = np.mean(dynrng[ix[i]:ix[i+1]-1])

    return Dout


############## Provide a command-line wrapper

# For SRI's wavreading code
import scipy.io.wavfile as wav
from scikits.audiolab import Sndfile
# For command line
import os
import sys

def readsph(filename):
    """ read in audio data from a sphere file.  Return d, sr """
    f = Sndfile(filename, 'r')
    data = f.read_frames(f.nframes, dtype=np.float32)
    sr = f.samplerate
    return data, sr

def readwav(filename):
    """ read in audio data from a wav file.  Return d, sr """
    # Read in wav file
    sr, wavd = wav.read(filename)
    # normalize short ints to floats of -1 / 1
    data = np.asfarray(wavd) / 32768.0  
    return data, sr

def audioread(filename, targetsr=None):
    """
    Read a soundfile of either WAV or SPH, based on filename
    returns d, sr
    """
    fileName, fileExtension = os.path.splitext(filename)
    if fileExtension == ".wav":
        data, sr = readwav(filename)
    elif fileExtension == ".sph":
        data, sr = readsph(filename)
    else:
        raise NameError( ("Cannot determine type of infile " +
                          filename) )
    # Maybe fix sample rate
    #if srate == 16000 and self.sbpca.srate == 8000:
    if targetsr != None and sr != targetsr:
        # Right now, only downsample by integer numbers
        decimfact = int(np.round(sr/targetsr))
        data = scipy.signal.decimate(np.r_[data[1:], 0], 
                                     decimfact, ftype='fir')
        # slight trim to ss.decimate to make its phase align 
        # to matlab's resample 
        # for case of resampling 16 kHz down to 8 kHz
        delay = 7
        data = np.r_[data[delay:], np.zeros(delay)]
        sr = sr/decimfact

    return data, sr

def main(argv):
    """ Main routine to apply from command line """
    if len(argv) != 2:
        raise NameError( ("Usage: " + argv[0] + 
                          " inputsound.{wav,sph}") )

    inwavfile = argv[1]

    data, sr = audioread(inwavfile, targetsr=8000)

    # Apply
    ftrs = ehist_dynrange(data, sr)

    # Write the data out
    print "%s %.1f %.1f %.1f %.1f %.1f" % (inwavfile, ftrs[0], ftrs[1], ftrs[2], ftrs[3], ftrs[4])

# Actually run main
if __name__ == "__main__":
    main(sys.argv)

