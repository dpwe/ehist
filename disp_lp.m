function [Xo, CVo, No] = disp_lp(name, pane, randomize)
% [envTypeCode, datafiles] = disp_lp(name, pane, randomize)
%   Display the energy histogram statistics for some utterance
%   set (language pack).  Usages:
%   disp_lp({'BP_101', 'BP_105', 'BP_106'})
%     - plot summary statistics of several entire language packs
%   disp_lp('BP_101') - an entire LP, broken down by envType
%   disp_lp({'BP_101_CK', 'BP_101_HOL'}) - summary stats of several conditions
%   disp_lp('BP_101_CK') - display sampling of utterances for a
%                     condition, spanning full range of WERs (if available)
%   disp_lp({'BABEL_OP1_206_65882_20121201_174526_outLine', ..})
%     - plot individual utterances
%   In each case, data is read from precalculated feature files
%   from calc_ehist_ftrs().
%   Return vector of envType assigments, and corresponding file
%   names.
%   <pane> set means to plot only into a single pane (1:8), and
%   don't divide across multiple conditions (where available).
%   <randomize> set shuffles the envType codes, just for a sanity check.
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; pane = 0; end
if nargin < 3; randomize = 0; end

maxnpanes = 8;

if iscell(name)
  % If we're given a list, it means to distribute them across panes
  Nn = length(name);
  % ensure existing panes are cleared
  clf;
  for i = 1:Nn
    pane = i;
    [XX,CC,NN] = disp_lp(name{i}, pane, randomize);
    if i == 1
      X{1} = XX{1}; CV{1} = CC{1}; N = NN;
    else
      X{1} = X{1} + XX{1}; CV{1} = CV{1} + CC{1}; N = N + NN;
    end
  end

  X{1} = X{1}/Nn;
  CV{1} = CV{1}/Nn;
  
  % Trick final plot, falls through
  nconds = 1;
  pane = Nn+1; % final pane
  titles{1} = 'ALL';
  
else
  % name is a singleton, have to figure out what it means
  % could be corpus 'OP1_206', corpus + condition 'OP1_206_CK', or
  % individual utterance 'BABEL_OP1_206_65882_20121201_174526_outLine'

  % is this an utterance?
  if length(name) > 6 && strcmp(name(1:6), 'BABEL_')
    % single utterance
    [X{1},CV{1},wer,cond] = ftrs_for_utt(name);
    nconds = 1;
    N = 1;
    titles{1} = [name, ' - ', cond];
  
  else
    underscores = find(name == '_');
    if length(underscores) ~= 1 && length(underscores) ~= 2
      disp(['Unable to parse name ', name]);
      return;
    end
    if length(underscores) == 1
      % Name is bare corpus - 'BP_101'
      corpus = name;
      conditions = [];
    else % length(underscores)==2, two underscores - 'BP_101_CK'
      corpus = name(1:(underscores(2)-1));
      conditions{1} = name(underscores(2)+1:end);
    end
 
    % Figure what utterances we have
    uttdir = babel_corpusdir(corpus);
    utts = dir(fullfile(uttdir, '*.mat'));
    % Convert dir output to list of uttids
    for i = 1:length(utts)
      % lower case and remove '.mat'
      dots = max(find(utts(i).name == '.'));
      %uttids{i} = lower(utts(i).name(1:(dots-1)));
      % No lower case .. ?
      uttids{i} = utts(i).name(1:(dots-1));
    end
    disp([corpus,': ',num2str(length(uttids)), ...
          ' utterance feature files found']);

    % Read the metadata for this corpus, including these utts
    [dem, asr, envTypeCode, wers] = read_ref_mats(corpus, uttids);
    
    envNames = dem.envTypes.abbrev;

    % Convert conditions to indices
    if length(conditions) == 0
      conditions = 1:length(envNames);
    else
      conditionNames = conditions;
      conditions = [];
      for i = 1:length(conditionNames);
        ix = strmatch(conditionNames{i}, envNames, 'exact');
        if length(ix) == 0
          disp(['unable to find condition ', conditionNames{i}]);
        else
          conditions(i) = ix;
        end
      end
    end
    
    if randomize
      % RANDOMIZE - CONTROL
      disp('randomizing condition assignments')
      envTypeCode = envTypeCode(randperm(length(envTypeCode)));
    end

    if pane ~= 0 && length(conditions) > 1
      % Collapse whole corpus to one condition if we only have a
      % single pane to display them in
      envTypeCode(:) = 1;
      conditions = [1];
      envNames{1} = '';
    end

    nconds = length(conditions);

    if pane == 0 && length(conditions) == 1
      % We're given a single condition to array across all panes
      % Subselect from utterances, sorted by WER
      uix = find(envTypeCode == conditions(1));
      % Sort by wer
      if length(wers)
        [vv,xx] = sort(wers(uix), 'descend');
        uix = uix(xx);
      end
      % stride down uix list
      uixix = 1 + round((length(uix)-1)*[0:maxnpanes-1]/(maxnpanes-1));
      for i = 1:maxnpanes
        N(i) = 1;
        uttix = uix(uixix(i));
        uttid = uttids{uttix};
        [X{i}, CV{i}, wer(i)] = ftrs_for_utt(uttid, dem, asr);
        titles{i} = [uttid, ' - ', envNames{envTypeCode(uix(uixix(i)))}];
      end
      nconds = maxnpanes;
    else
      % Figure out the contents for each pane by pooling the
      % utterances matching each condition
      for i = 1:nconds
        uix = find(envTypeCode == conditions(i));
        N(i) = length(uix);
        [X{i}, CV{i}, wer(i)] = ftrs_for_utt(uttids(uix), dem, asr);
        titles{i} = [name, ' - ', envNames{i}];
      end
    end
  end
end

% At this point, X{i} and CV{i} are the (possibly pooled) energy
% histograms and covariance matrices to display in each pane i,
% with associated utterance counts N(i) and overall WER wer(i).  

if pane == 0; 
  % first call of multi-pane plot
  pane = 1; 
  % ensure existing panes are cleared
  clf;
end

for i = 1:length(N)
  if N(i) > 0
    if N(i) > 1
      subtitle = ['N=', num2str(N(i))];
    elseif wer(i) > 0
      subtitle = ['WER=', sprintf('%.1f', 100*wer(i))];
    else
      subtitle = '';
    end
    disp_pane(X{i}, CV{i}, pane, maxnpanes, titles{i}, subtitle);
  end
  pane = pane + 1;
end

% Assign outputs? (used to calculate summary pane in multi-pane displays)
if nargout > 0
  Xo = X;
  CVo = CV;
  No = N;
end
