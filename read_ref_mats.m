function [dem, asr, condcodes, wers] = read_ref_mats(corpus, utts)
% [dem, asr, condcodes, wers] = read_ref_mats(corpus, utts)
%   Read the reference materials (demographics, ASR results if
%   available) for named corpus.
% 2014-01-09 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; utts = []; end

% Find which conditions they have
demdir = 'reference_materials';
demographicsfile = fullfile(demdir, [corpus, '-demographics.tsv']);
asrresultsfile   = fullfile(demdir, [corpus, '-dev.ctm.sys']);
dem = demographics_read(demographicsfile);
%disp([num2str(length(dem.outputFn)), ' demographics records read']);
asr = asrresults_read(asrresultsfile, dem);
%disp([num2str(length(asr.uttID)), ' ASR results read']);

nutts = length(utts);

if nutts > 0
  wers = zeros(1, nutts);
  condcodes = zeros(1, nutts);
  dixes = zeros(1, nutts);
  for i = 1:nutts
    %% remove .mat
    %uttID = lower(utts(i).name(1:end-4));
    uttID = utts{i};
    
    dix = strmatch(lower(uttID), lower(dem.outputFn), 'exact');
    % record which condition this utt belongs to
    if length(dix)
      condcodes(i) = dem.envTypes.code(dix);
      dixes(i) = dix;
    end
    aix = strmatch(uttID, asr.uttID, 'exact');
    if length(aix)
      % WER for this utt
      wers(i) = asr.err(aix)/asr.chr(aix);
    end
  end
  % subselect the relevant parts of dem
  fldnames = fieldnames(dem);
  for fldix = 1:length(fldnames)
    field = fldnames{fldix};
    val = getfield(dem, field);
    if size(val, 2) > 1
      dem = setfield(dem, field, val(dixes));
    end
  end
else
  condcodes = [];
  wers = [];
end
