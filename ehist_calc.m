function [X,DNM,CV,F] = ehist_calc(N, domask, offset, range, melbins)
% [X,DNM,CV,F] = ehist_calc(N, domask, offset, range, melbins)
%    X return a 2D histogram of levels/dB vs. freq for an entire
%    conversation side. 
%    offset is a dB offset added to every sgram value (dflt 70).
%    range is the largest dB value (scalar), or the range of dB
%    values to use (pair), or a set of dB bin centers (vector).
%    DNM returns the parent spectrogram, with noise gating if selected
%    Originally called uttftr.m .
%    CV returns a covariance matrix.
%    F returns the frequencies for each bin.
%    If melbins is greater than zero, map into mel bins instead of
%    sgram bins.
% 2013-05-25 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; domask  = 1;   end
if nargin < 3; offset  = 80;  end
if nargin < 4; range   = 120; end
if nargin < 5; melbins = 0;   end

targetsr = 8000;
forcemono = 1;

[dn,sr] = audioread(N,targetsr,forcemono);

%fftsize = 4096;
%
%DN = abs(specgram(dn,fftsize,sr));
%
%% average and smooth
%
%mD = mean(DN');
%cD = cov(DN');
%
%smwin = 15;
%
%smD = conv(mD, hann(smwin)/sum(hann(smwin)));
%smD = smD(floor(smwin/2)+[1:(fftsize/2+1)]);
%
%X = smD;

% from snreval.m

dlen = length(dn);
DISP = 0;

%dbmin = 0.0;
%dbmax = 120.0;
dbbinwidth = 1.0;
if length(range) == 1
  range = [0 range];
end
if length(range) > 2
  dbbins = range;
else
  dbbins = range(1):dbbinwidth:range(2);
end

% Choose fft size; 1024 for sr = 16000
nfft = 2^round(log(1024 * sr/16000)/log(2));
nhop = nfft/2;
nsgframes = 1+ floor((dlen-nfft)/nhop);
fr = sr/nhop;
% Signal spectrogram
SG = abs(specgram(dn,nfft));
% lose the nyquist bin as it messes up the histogram
SG = SG(1:size(SG,1)-1,:);

if melbins > 0
  width = 2.0;  % double-width bins
  [wts, F] = fft2melmx(nfft, sr, melbins, width);
  DN = 20*log10(wts(:,1:(nfft/2))*SG);
else
  DN = 20*log10(SG);
  F = [0:(nfft/2)]/nfft*sr;
end

% collapse to 0..80 dB range (from top)
%maxdn = max(DN(:));
%DN = max(0,dbmax+DN-max(DN(:)));
DN = DN + offset;

% Build a sad mask from a dB spectrogram
percentl = 0.1;
marg     = 10.0;
collar   = 5; 
SM = sadmask(DN, percentl, marg, collar);
if ~domask
  SM = 1+0*SM;
end
% Keep only the active parts
DNM = DN(:, find(SM));
% now index into full-length DN

% Calculate subband energy histogram in each row, in 1 dB bins
X = zeros(length(dbbins), size(DNM,1));
for i = 1:size(DN,1);
  X(:,i) = hist(DNM(i,:),dbbins)';
end

% Covariance matrix
%CV = 10*log10(cov(exp(DNM'/8.68588963806504)));
CV = (cov(exp(DNM'/8.68588963806504)));


% then:
% >> bp101root = ...
% '/u/drspeech/data/swordfish/corpora/BABEL_BP_101/conversational';
% >> bp101dir = fullfile(bp101root, 'training/audio');
% >> [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12] = ...
%      textread(fullfile(bp101root,'reference_materials/demographics.tsv'), ...
%      '%s%s%s%s%s%s%s%s%s%s%s%s','whitespace','\t');
% >> % strip header line
% >> wavs = m1(2:end);
% >> conds = m9(2:end);
% >> nwav = length(wavs);
% >> XX = zeros(nwav,121,257);
% >> for i = 1:nwav; disp(wavs{i}); XX(i,:,:) = uttftr(fullfile(bp101dir,wavs{i})); end
% >> [C,IA,IC] = unique(conds);  % IC is condition index per wav
% >> GM = zeros(length(C),size(XX,2),size(XX,3));
% >> for i = 1:nwav; GM(IC(i),:,:) = GM(IC(i),:,:) + XX(i,:,:); end
% >> sr = 8000;
% >> ff = [0:size(XX,3)-1]/size(XX,3)*sr/2;
% >> for i = 1:length(C); subplot(3,2,i); imgsc(ff,0:120,squeeze(GM(i,:,:))); title([conds{IA(i)},' - ',num2str(sum(IC==i))], 'interpreter', 'none'); caxis([0 2*max(max(squeeze(GM(i,30:100,100:200))))]); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
