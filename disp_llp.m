function [envTypeCode, datafiles, wer, X, CV, N] = disp_llp(name, noconds, pane, randomize)
% [envTypeCode, datafiles] = disp_llp(name, pane, noconds, randomize)
%   Display the freq statistics for an entire LLP under <name>
%   (e.g 'BP_101')
%   by reading in all the precalculated utterance feature files
%   (from calc_ftrs) and accumulating them according to the
%   metadata (from condMap).
%   Return vector of envType assigments, and corresponding file
%   names.
%   <nocond> set means to suppress the final summary pane.
%   <pane> set means to plot only into a single pane (1:8).
%   <randomize> set shuffles the envType codes, just for a sanity check.
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; pane = 0; end
if nargin < 3; noconds = 0; end
if nargin < 4; randomize = 0; end

if iscell(name)
  % list of dirs - plot them all
  Nn = length(name);
  for i = 1:Nn
    pane = i;
    [ee,dd,ww,XX,CC,NN] = disp_llp(name{i}, noconds, pane, randomize);
    if i == 1
      X{1} = XX{1}; CV{1} = CC{1}; N = NN;
    else
      X{1} = X{1} + XX{1}; CV{1} = CV{1} + CC{1}; N = N + NN;
    end
  end

  % Trick final plot, falls through
  nconds = 1;
  noconds = length(name)+1; % final pane
  name = 'Corpora';
  envNames{1} = 'ALL';
  
else
  % normal call - single corpus, break down by conditions

  corpus = name;
  
  % Figure what utterances we have
  uttdir = ['corpora/BABEL_', corpus, ...
            '_LLP/conversational/dev/audio'];
  utts = dir(fullfile(uttdir, '*.mat'));
  disp([num2str(length(utts)), ' utterance feature files found']);

  % Find which conditions they have
  demdir = 'wer_by_condition';
  demographicsfile = fullfile(demdir, [corpus, '-demographics.tsv']);
  asrresultsfile   = fullfile(demdir, [corpus, '-dev.ctm.sys']);
  dem = demographics_read(demographicsfile);
  disp([num2str(length(dem.outputFn)), ' demographics records read']);
  asr = asrresults_read(asrresultsfile, dem);
  disp([num2str(length(asr.uttID)), ' ASR results read']);

  envNames = dem.envTypes.names;

  % create env index for utterances
  for i = 1:length(utts)
    % remove .mat
    uttID = lower(utts(i).name(1:end-4));
    ix = strmatch(uttID, dem.outputFn);
    % record which condition this utt belongs to
    envTypeCode(i) = dem.envTypes.code(ix);
    aix = strmatch(uttID, asr.uttID);
    if length(aix)
      % WER for this utt
      wer(i) = asr.err(aix)/asr.chr(aix);
      cond = asr.envTypeNames{asr.envTypeCode(aix)};
    else
      wer(i) = 0;
      cond = '';
    end
    % also set up simple cell array of full file names
    datafiles{i} = fullfile(uttdir, utts(i).name);
  end

  if randomize
    % RANDOMIZE - CONTROL
    disp('randomizing condition assignments')
    envTypeCode = envTypeCode(randperm(length(envTypeCode)));
  end

  if noconds
    % Collapse whole corpus to one condition
    envTypeCode(:) = 1;
    envNames{1} = 'ALL';
  end

  warned = 0;

  nconds = max(envTypeCode);

  for i = 1:nconds
    uix = find(envTypeCode == i);
    N(i) = length(uix);
    if N(i) == 0
      X{i} = [];
      CV{i} = [];
    else
      for j = 1:N(i)
        D = load(datafiles{uix(j)});
        if max(abs(D.CV(:))) < 100
          if ~warned; disp('Fixing log CV'); warned = 1; end
          D.CV = real(exp(log(10)/10*D.CV));
        end
        if j == 1
          X{i} = D.X;
          CV{i} = D.CV;
        else
          X{i} = X{i} + D.X;
          CV{i} = CV{i} + D.CV;
        end
      end
    end
  end
end

pane = 0;
for i = 1:nconds
    
  if N(i) > 0
    if noconds == 0
      pane = pane + 1;
    else
      pane = noconds;
    end
    subplot(4,4,2*pane-1)
    tit = [name, ' - ', envNames{i}];
    ehist_disp(X{i}/N(i), [], 0, tit, (pane>4), mod(pane,2) );
    % Add count of utterances
    textypos = 0;
    textxpos = 1;
    text(textxpos, textypos, num2str(N(i)), 'Color', [0.9 0.9 0.9], ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

    % also cov
    subplot(4,4,2*pane);
    ehistcv_disp(CV{i}/N(i), '', pane>4);

    end
  end
end
