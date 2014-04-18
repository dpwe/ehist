#!/usr/bin/env python
#
# ehist_enhance.py
#
# Enhance an utterance by equalizing its subband energy histograms.
#
#2014-04-17 Dan Ellis dpwe@ee.columbia.edu

import numpy as np
# For filter
import scipy.signal
# mel mapping
import fft2melmx

############### stft analysis/synthesis ##################

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

def ola(X, hop):
    """ Overlap-add rows of X into a single vector, stepping by hop """
    nw, W = np.shape(X)
    # How long should X be in the end?
    # one complete window, plus another hop's worth of points for each
    # additional window
    lenx = W + (nw-1)*hop
    Y = np.zeros(lenx)

    for i in range(nw):
        Y[i*hop : i*hop+W] += X[i,]

    return Y

def stft(signal, nfft, windowlen, hop):
    """ calculate the short-time fourier transform """
    frames = frame(signal, windowlen, hop)

    # apply frame window to each frame
    window = np.hanning(windowlen)
    wframes = frames * window

    return np.fft.rfft(wframes, int(nfft))

def stftm(signal, nfft, windowlen, hop):
    """ calculate the short-time fourier transform magnitude """
    return np.abs(stft(signal, nfft, windowlen, hop))

def istft(ffts, windowlen, hop):
    """ undo the effect of stft.
        Set windowlen to 0 to defeat re-windowing """
    frames = np.fft.irfft(ffts)
    if windowlen > 0:
        window = np.hanning(windowlen)
        # .. centered extension/trimming here
        frames *= window
    return ola(frames, hop)

################ piecewise-linear-mapping routines ##############

def map_vals(X, vmap):
    """ X is a vector, map is two rows indicating input, output pairs """
    if len(np.shape(X)) > 1:
        # If X is a matrix, run separately on each row
        return np.array([map_vals(XX, vmap) for XX in X])
    else:
        # Just run on a vector
        # figure slopes
        gap = np.diff(vmap)
        # repeat final gap on both src and dst rows
        gap = np.c_[gap, gap[:,-1][:, np.newaxis]]
        # Don't allow zero (or negative?) gaps
        gap[np.nonzero(gap <= 0)] = np.finfo(float).eps
        # do mapping
        Xix = np.maximum(0, 
                         np.sum(np.greater.outer(X, vmap[0]), axis=1)
                         - 1 )
        Xdelta = (X - vmap[0, Xix])/gap[0, Xix]
        return vmap[1, Xix] + Xdelta*gap[1, Xix]

def inv_map(vmap):
    """ construct the inverse of a monotonic map s.t. 
        map_vals(map_vals(vals, vmap), inv_map(vmap)) == vals """
    # you simply interchange the input points and the output points
    return vmap[::-1]

def compose_maps(map1, map2):
    """ return a single [x,y] map that represent y = map1(map2(x)) """
    # Break points will be all the edges in map2, and all the edges in map1 
    # when projected through the inverse of map2
    mapped_edges = map_vals(map1[0], inv_map(map2))
    all_edges = np.unique(np.r_[map2[0], mapped_edges])
    all_vals = map_vals(map_vals(all_edges, map2), map1)
    return np.array([all_edges, all_vals])

############# mapping of histograms ###################

def make_hist_map(hist_in, hist_out, edges):
    """ return a mapping that will map the values from hist_in
        to match those in hist_out """
    # for each row of the histograms
    #   cumsum & normalize
    #   figure where in one the other occurs
    cdf_in = np.cumsum(hist_in)/np.sum(hist_in)
    cdf_out = np.cumsum(hist_out)/np.sum(hist_out)
    cdf_in_map = np.array([edges, cdf_in])
    cdf_out_map = np.array([edges, cdf_out])
    # The histogram map is
    # x_mapped = cdf_out^{-1}(cdf_in(x_in))
    return compose_maps(inv_map(cdf_out_map), cdf_in_map)

def make_hist_maps(hists_in, hists_out, edges):
    """ make a list of maps, one for each matching row of hists_in and out """
    return [make_hist_map(hist_in, hist_out, edges) 
            for hist_in, hist_out in zip(hists_in, hists_out)]

def apply_hist_maps(X, histmap):
    """ map values in X according to the histogram map
        histmap is a set of pairs defining piecewise constant maps
        <inval, outval>
        histmap may be a list of maps, each applied to a row of X 
    """
    return np.array([map_vals(row, vmap) for row, vmap in zip(X, histmap)])

################# histograms and percentiles ###################

def percentile(data, pcntl):
    """
    v = percentile(d,n)
    Return for each col of v the n'th percentile where 0<n<1
    2004-10-04 dpwe@ee.columbia.edu
    """
    nr = np.size(data, axis=0)

    x = np.sort(data, axis=0)
    return x[int(np.floor(pcntl*nr)),]

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
    # discard top bin
    if np.sum(counts[:,-1]) > 0:
        print "histc: warning: some values above highest edge"
    return counts[:,:-1]

############## deciBels ##################

def dB(X):
    """ map linear amplitude to dB """
    return 20*np.log10(X)

def idB(Y):
    """ map dB to linear amplitude """
    return np.exp(Y/8.685889638065035)

############### mel_hist ###################

def mel_hist(d, sr, nfft=256, win=None, hop=None, nfilts=40, mindb=-100., maxdb=0., dbbin=1.0, edges=None):
    """ Calculate the mel-freq energy histogram for some audio """
    # Base stft
    if win == None:
        win = nfft
    if hop == None:
        hop = win/4
    D = stft(d, nfft, win, hop).T
    # Map to mel
    melmx, freqs = fft2melmx.fft2melmx(nfft, sr, nfilts)
    DmeldB = dB(np.dot(melmx, np.abs(D)))
    # Get histogram in mel bins
    if edges == None:
        edges = np.arange(mindb, maxdb, dbbin)
    melHist = histc(DmeldB, edges)
    # discard lowest bin of histogram (underflow)
    if np.sum(melHist[:,0]) > 0:
        print "mel_hist: Warning: discarding %f of undeflow" % np.sum(melHist[:,0])
        melHist[:,0] = 0

    return melHist, edges, D, DmeldB, melmx, freqs

########### actually modify signal to map mel histograms ########
from scipy.ndimage.filters import median_filter

def ehist_equalize_melhist(d, sr, refMelHist, edges):
    """ Modify a signal in the Mel domain 
        by equalizing the Mel-subband histograms to match
        the passed-in ones """
    # Calculate the (Mel) spectrograms, and histogram, and axes
    melHist, edges, D, DmeldB, melmx, freqs = mel_hist(d, sr, edges=edges)
    # Build mapping & modify mel spectrogram
    histmaps = make_hist_maps(melHist, refMelHist, edges)
    # for some reason, extrapolating madly below bottom edge - clip it
    DmeldBmapped = np.maximum(edges[0], 
                              np.minimum(edges[-1], 
                                         apply_hist_maps(DmeldB, histmaps)))
    # Reconstruct audio based on mapped envelope
    # We map both original and modified Mel envelopes to FFT domain
    # then scale original STFT magnitudes by their ratio
    DmelInFFT = np.dot(melmx.T, idB(DmeldB))
    DmappedInFFT = np.dot(melmx.T, idB(DmeldBmapped))
    # Zero values in denominator will match to zeros in numerator, 
    # so it's OK to drop them
    Dmask = DmappedInFFT / (DmelInFFT + (DmelInFFT==0))
    # Median filter to remove short blips in gain
    medfiltwin = 7
    DmaskF = median_filter(Dmask, size=(1, medfiltwin))
    # Now scale each FFT val by their ratio
    Dmod = D * DmaskF
    # and resynthesize
    nfft = 2*(np.size(D, axis=0)-1)
    win = nfft
    hop = win/4
    dout = istft(Dmod.T, win, hop)
    return dout

########### main function ##############

def ehist_enhance(infile, inref, outfile):
    """
    Enhance the input file by matching its mel-subband energy histograms 
    to those of the reference file.  Write the output file
    """
    # build the reference histogram
    dref, srref = audioread(inref)
    refHist, edges, Dr, Drmel, melmx, freqs = mel_hist(dref, srref)
    # read in the data to modify
    d, sr = audioread(infile)
    # modify
    dmod = ehist_equalize_melhist(d, sr, refHist, edges)
    # save
    audiowrite(dmod, sr, outfile)

############## Provide a command-line wrapper

# For SRI's wavreading code
import scipy.io.wavfile as wav
from scikits.audiolab import Sndfile, Format
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

def audiowrite(data, sr, filename):
    """
    Write audio data to a file.  Infer type from name.
    """
    stem, ext = os.path.splitext(filename)
    format = Format(ext[1:])
    if len(np.shape(data)) > 1:
        nchans = np.size(data, axis = 0)
    else:
        nchans = 1
    f = Sndfile(filename, 'w', format, nchans, sr)
    if np.max(np.abs(data)) >= 0.999:
        clippos = data >= 0.999
        data[clippos] = 0.999
        clipneg = data <= -0.999
        data[clipneg] = -0.999
        print "audiowrite: WARNING: %d samples clipped" % np.sum(np.r_[clippos, clipneg])
    f.write_frames(data)
    f.close()

def main(argv):
    """ Main routine to apply from command line """
    if len(argv) != 4:
        raise NameError( ("Usage: " + argv[0] + 
                          " inputsound.wav refsound.wav outsound.wav") )
    insound = argv[1]
    refsound = argv[2]
    outsound = argv[3]
    ehist_enhance(insound, refsound, outsound)


# Actually run main
if __name__ == "__main__":
    main(sys.argv)

