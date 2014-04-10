function [asr, dem] = analyze_conditions(corpus, randomize, bars)
% [asr, dem] = analyze_conditions(corpus, randomize, bars)
%   For a corpus (e.g. 'BP_101') read an asrresults and a
%   demographics file, and calculate the results per condition.
%   <randomize> as 1 shuffles the condition labels as a control.
%   <bars> selects stacked bars instead of boxplot.
% 2014-01-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; randomize = 0; end
if nargin < 3; bars = 0; end

if iscell(corpus)
  % multiple corpora - plot them all
  for i = 1:length(corpus); 
    subplot(3,2,i); 
    a{i} = analyze_conditions(corpus{i}, randomize, bars); 
  end
  % add summary pane
  asr = a{1};
  for i = 2:length(corpus); 
    asr.chr = [asr.chr, a{i}.chr]; 
    asr.snt = [asr.snt, a{i}.snt]; 
    asr.corr = [asr.corr, a{i}.corr]; 
    asr.sub = [asr.sub, a{i}.sub]; 
    asr.ins = [asr.ins, a{i}.ins]; 
    asr.del = [asr.del, a{i}.del]; 
    asr.err = [asr.err, a{i}.err]; 
    asr.serr = [asr.serr, a{i}.serr];
    asr.typeCode = [asr.typeCode, a{i}.typeCode]; 
  end
  subplot(326)
  corpus = 'Sum';  % name for final plot
  % .. then fall through to boxplot command
  addlegend = 1;
  
else
  % just plotting a single corpus
  
  [dem, asr] = read_ref_mats(corpus);
  % Which condition do we factor across?
  % envType
  asr.typeCode  = asr.envTypeCode;
  asr.typeNames = asr.envTypeNames;
  % gender
%  asr.typeCode  = asr.gender;
%  asr.typeNames = asr.genderNames;
  % dialect
%  asr.typeCode  = asr.dialectCode;
%  asr.typeNames = asr.dialectNames;
  
  if randomize
    % RANDOMIZE - CONTROL
    disp('randomizing type assignments')
    asr.typeCode = asr.typeCode(randperm(length(asr.typeCode)));
  end
  
  addlegend = 0;

end  

if bars
  plot_conditions_bar(asr, corpus, addlegend);
else
  plot_conditions_box(asr, corpus, addlegend);
end

% To generate plot:
% corpora = {'BP_101', 'BP_104', 'BP_105', 'BP_106', 'BP_107'};
% analyze_conditions(corpora);
% print -depsc bp_1xx_err-by-condition.eps
