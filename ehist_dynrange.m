function [D,Nout] = ehist_dynrange(N)
% [D,Nout] = ehist_dynrange(N)
%   Calculate an energy histogram for the specified soundfile
%   and report the dB difference between 5th and 95th percentiles
%   for each of the 5 octaves 
%   125-250 Hz, 250-500 Hz, 500-1 kHz, 1-2 kHz, 2-4 kHz.
%   If N is a cell array, apply to each, return N x 5 matrix.
%   If N names a directory, apply to all *.sph in that directory.
%   Nout lists the actual files processed.
% 2014-04-09 Dan Ellis dpwe@ee.columbia.edu

% decoration...

if iscell(N)
  % Apply to each element of a cell array
  for i = 1:length(N)
    [D(i,:), Nout{i}] = ehist_dynrange(N{i});
  end
elseif isdir(N)
  % Passed a directory name - apply to all files in that directory
  uttdir = N;
  utts = dir(fullfile(uttdir, '*.sph'));
  N = cell(length(utts));
  for i = 1:length(utts)
    N{i} = fullfile(uttdir, utts(i).name);
  end
  [D, Nout] = ehist_dynrange(N);
else
  
  % Assume passed an actual wav file name

  domask = 1;  % Apply VAD masking (why?)
  offset = 90; % shift up by 90 dB
  range  = 120; % report 120 x 1 dB steps
  melbins = 0;  % report in FFT bins, don't map to Mel

  [X,DNM,CV,F] = ehist_calc(N, domask, offset, range, melbins);

  % Calculate difference of 5th and 95th percentiles
  pctls = [0.05 0.95];
  dynrng = diff(histpercentile(X, pctls));

  % Find bin edges
  fmin = 125;
  fmax = 4000;
  fthis = fmin;
  for i = 1:round(1+log(fmax/fmin)/log(2))
    ix(i) = max(find(F <= fthis));
    fthis = 2*fthis;
  end

  % Calculate average dyn ranges in those fbin ranges
  for i = 1:length(ix) - 1
    Dout(i) = mean(dynrng(ix(i):ix(i+1)-1));
  end

  if nargout < 1
    % just report result
    for i = 1:length(Dout)
      disp(['Dyn range ', sprintf('%.0f - %.0f Hz = %.1f dB', ...
                                  F(ix(i)), F(ix(i+1)), Dout(i))]);
    end
  else
    D = Dout;
    Nout = N;
  end

end
