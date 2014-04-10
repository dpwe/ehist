function [X,CV,wer] = disp_utts(utts, pane)
% [X, CV, wer] = disp_utts(utts, pane)
%   Display the freq statistics for a set of individual tracks.
%   Include ASR performance where available.
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; pane = 0; end

if iscell(utts)
  % list of utts - plot them all
  N = length(utts);
  for i = 1:N
    subpane = i;
    [XX,CC,ww] = disp_utts(utts{i}, subpane);
    if i == 1
      X = XX; CV = CC; wer = ww;
    else
      X = X + XX; CV = CV + CC; wer = wer + ww;
    end
  end

  
  X = X/N;
  CV =CV/N;
  wer = wer/N;

  if pane < 0
    % *don't* plot the summary
    return
  end

  % now fall through & plot sum
  pane = 6;
  utt = 'Sum';
  cond = '';
  
else
  % Single utterance
  utt = utts

  [p,n,e] = fileparts(utt);

  utt = [n,e];
  corpus = utt([7:13]);
  if corpus(end) == '_'; corpus = corpus(1:end-1); end
    
  if length(p) == 0
    % Default path to LLP dev utts
    uttdir = ['corpora/BABEL_',corpus, ...
              '_LLP/conversational/dev/audio'];
  else
    uttdir = p;
  end
  
    
  % Find which conditions they have
  demdir = 'wer_by_condition';
  demographicsfile = fullfile(demdir, [corpus, '-demographics.tsv'])
  asrresultsfile   = fullfile(demdir, [corpus, '-dev.ctm.sys']);
  dem = demographics_read(demographicsfile);
  disp([num2str(length(dem.outputFn)), ' demographics records read']);
  asr = asrresults_read(asrresultsfile, dem);
  disp([num2str(length(asr.uttID)), ' ASR results read']);

  envNames = dem.envTypes.names;

  % remove .mat
  uttID = lower(utt);
  dix = strmatch(uttID, dem.outputFn);
  % record which condition this utt belongs to
  cond = '';
  if length(dix)
    envTypeCode = dem.envTypes.code(dix);
    cond = dem.envTypes.abbrev{envTypeCode};
  end
  % Look it up in the ASR results
  aix = strmatch(uttID, asr.uttID);
  wer = 0;
  if length(aix)
    % WER for this utt
    wer = asr.err(aix)/asr.chr(aix);
  end
  % also set up simple cell array of full file names
  datafile = fullfile(uttdir, [utt, '.mat']);

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

if pane
  subplot(4,4,2*pane-1);
end
tit = [utt, ' - ', cond];
ehist_disp(X, [], 0, tit, pane>6, mod(pane,2));
if wer > 0
  % Add count of utterances
  textypos = 0;
  textxpos = 1;
  text(textxpos, textypos, ['WER=',sprintf('%.1f', 100*wer)], ...
       'Color', [1 1 1], ...
       'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

% also cov
if pane
  subplot(4,4,2*pane);
  ehistcv_disp(CV, '', pane>6);
end
