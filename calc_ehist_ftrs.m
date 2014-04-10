function calc_ehist_ftrs(N,dstdir)
% calc_ehist_ftrs(N,dstdir)
%   Calculate the utterance energy histogram features for file N or
%   list of files. 
%   Save to dstdir/N.mat the histogram and the covariance matrix. 
% 2014-01-05 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; dstdir = 'ftrsout'; end

if iscell(N)
  % Run over each item in a cell array
  nitems = length(N);
  for i = 1:length(N)
    calc_ehist_ftrs(N{i}, dstdir);
  end

elseif isdir(N)
  % Run on all the sph files in that directory
  wavdir = N;
  fdir = dir(fullfile(wavdir, '*.sph'));
  for i = 1:length(fdir) 
    fl{i} = fullfile(wavdir, fdir(i).name); 
  end
  calc_ehist_ftrs(fl, dstdir);
  
else
  % Just a file name to process
  [p,n,e] = fileparts(N);
  outdir = fullfile(dstdir, p);
  mymkdir(outdir);
  outpath = fullfile(outdir, [n,'.mat']);

  domask = 1;
  offset = 90;
  range = 120;
  melbins = 80;
  [X,DNM,CV,F] = ehist_calc(N, domask, offset, range, melbins);

  save(outpath, 'X', 'CV');
  disp(['wrote ', outpath]);

end

% 2014-01-05 Creating feature files for all LLP training utts on login.icsi.berkeley.edu
%
% llps = {'BP_101','BP_104','BP_105','BP_106', 'BP_107', 'OP1_102', 'OP1_103', 'OP1_201', 'OP1_203', 'OP1_206'}
% for i = 1:length(llps); llpp{i} = ['/u/drspeech/data/swordfish/corpora/BABEL_', llps{i},'_LLP/conversational/training/audio']; end
% tic; calc_ehist_ftrs(llpp); toc

