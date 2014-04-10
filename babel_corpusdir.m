function d = babel_corpusdir(corpus)
% d = babel_corpusdir(corpus)
%   Return the utterance directory name for a particular Babel corpus.
%   e.g. babel_corpusdir('BP_101') 
%   -> 'corpora/BABEL_BP_101_LLP/conversational/dev/audio'
% 2014-01-09 Dan Ellis dpwe@ee.columbia.edu

d = ['corpora/BABEL_', corpus, '_LLP/conversational/dev/audio'];
