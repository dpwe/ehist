function [X,CV,wer,cond, condix] = ftrs_for_utt(utts, dem, asr)
% [X, CV, wer, cond, condix] = ftrs_for_utt(utts, getwer)
%   Obtain the info to display for a single utt, or over a list.
%   Include ASR performance where available if getwer is set
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; dem = []; end
if nargin < 3; asr = []; end

cond = '';
condix = [];

if iscell(utts)
  % list of utts - plot them all
  N = length(utts);
  for i = 1:N
    [XX,CC,ww] = ftrs_for_utt(utts{i}, dem, asr);
    if i == 1
      X = XX; CV = CC; wer = ww;
    else
      X = X + XX; CV = CV + CC; wer = wer + ww;
    end
  end
  
  X = X/N;
  CV =CV/N;
  wer = wer/N;
  
else
  % Single utterance
  utt = utts;

  [p,n,e] = fileparts(utt);
  
  utt = [n,e];
  underscores = find(utt=='_');
  % Extract corpus name from utterance ID - between 1st and 3rd
  % underscores 
  corpus = utt( (underscores(1)+1):(underscores(3)-1) );
    
  if length(p) == 0
    % Default path to LLP dev utts
    uttdir = babel_corpusdir(corpus);
  else
    % Path was included as part of uttid
    uttdir = p;
  end
  
  % Find which conditions they have
  if length(dem) == 0
    [dem, asr] = read_ref_mats(corpus);
  end
  
  wer = 0;
  cond = '';
  condix = 0;
  if length(dem) > 0
    envNames = dem.envTypes.names;

    uttID = lower(utt);
    dix = strmatch(uttID, lower(dem.outputFn));
    % record which condition this utt belongs to
    if length(dix)
      condix = dem.envTypes.code(dix);
      cond = dem.envTypes.abbrev{condix};
    end
    % Look it up in the ASR results
    aix = strmatch(uttID, asr.uttID);
    if length(aix)
      % WER for this utt
      wer = asr.err(aix)/asr.chr(aix);
    end
    % also set up simple cell array of full file names
    datafile = fullfile(uttdir, [utt, '.mat']);
  end
  
  persistent warned;
  if ~exist('warned')
    warned = 0;
  end
  
  D = load(datafile);
  if max(abs(D.CV(:))) < 100
    if ~warned; disp('Fixing log CV'); warned = 1; end
    D.CV = real(exp(log(10)/10*D.CV));
  end
  X = D.X;
  CV = D.CV;

end
