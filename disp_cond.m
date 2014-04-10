function disp_cond(corpus, condition)
% disp_cond(corpus, condition)
%   Show some example utterances for a particular condition.
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

[envTypeCodes, datafiles, wers] = disp_llp(corpus);

ntoshow = 8;

% Find all relevant utterances
uix = find(envTypeCodes == condition);

% Sort by wer
[vv,xx] = sort(wers(uix), 'descend');
uix = uix(xx);

% Make uttID list
for i = 1:length(datafiles); 
  [p,utts{i},e] = fileparts(datafiles{i}); 
end

% display
disp_utts(utts(uix(1+round((length(uix)-1)*[0:ntoshow-1]/(ntoshow-1)))), ...
          -1); % no summary pane
