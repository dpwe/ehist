function [H,Bout] = colhist(X,Bin)
% [H,Bout] = colhist(X,Bin)
%    Form histograms of the data in each column of array X, using
%    the bins defined in Bin.  Bout gives the actual bin centers
%    used, and each column of H is the corresponding histogram of a
%    column of X.
%    Scalar Bin gives the number of bins to use.
% 2014-01-04 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; Bin = 10; end

if length(Bin) == 1
  % equally divide full range of data
  nbins = Bin;
  Xmin = min(X(find(~isnan(X))));
  Xmax = min(X(find(~isnan(X))));
  Xstep = (Xmax - Xmin)/(nbins-1);
  Bin = Xmin:Xstep:Xmax;
end

% Calculate histogram in each column
H = zeros(length(Bin), size(X,2));
for i = 1:size(H,2);
  [H(:,i), Bout] = hist(X(:,i), Bin);
end
