function M = sadmask(dBD, pcntl, marg, collar)
% M = sadmask(dBD, pcntl, marg, collar)
% Build a sad mask from a dB spectrogram
% pcntl is the point on the CDF to choose the low
% threshold (default 0.1 = 10th percentile).
% margin is the dB margin above the low threshold to count as the
% speech active threshold (dflt 10.0).
% collar is the number of frames to spread the speech active
% regions in both directions, and also median filter window 
% (dflt 10). 
% 2013-05-27, 2014-01-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; pcntl = 0.1; end
if nargin < 3; marg = 10.0; end
if nargin < 4; collar = 10; end

% total energy is taken as mean across all dB histo bins (!)
mdBD = mean(dBD);

thresh = percentile(mdBD, pcntl);
M1 = mdBD > (thresh + marg);

% grow it out this far in both directions
%collar = 10;  % with 16 ms hop, this is +/- 160ms
M = M1;
for i = 1:collar
  % OR together versions of the mask shifted both ways
  M = M | [zeros(1,i),M1(1:end-i)] | [M1((i+1):end),zeros(1, i)];
end

% Median filter too
%M = medianf(M,collar);
M = medfilt1(double(M),collar);

% but in all cases exclude frames with supernegative bins (digital
% silence)
M = M & (median(dBD)>0);
