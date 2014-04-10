function P = histpercentile(D,T)
% P = histpercentile(D,T)
%  Columns of D are histograms.  Return the 100*T'th percentile
%  (default 0.5) for each column in P.
%  If T is a vector, return multiple rows, one for each percentile.
% 2013-06-15 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; T = 0.5; end

nT = length(T);
P = zeros(nT, size(D,2));

oldway = 0;

if oldway
  
  % Calculate cumsums once per threshold
  for tx = 1:nT
    for i = 1:size(D,2)
      P(tx,i) = max(find(cumsum(D(:,i)) <= T(tx)*sum(D(:,i)))); 
    end
  end
  
else
  
  % Just one cumsum
  cs = cumsum(D)./repmat(sum(D),size(D,1),1);
  for tx = 1:nT
    P(tx,:) = sum(cs <= T(tx));
  end

end

% oldway:
%>> tic; PP = histpercentile(XC,[.05:.05:.95]); toc
%Elapsed time is 0.037798 seconds.
% newway:
%>> tic; PN = histpercentile(XC,[.05:.05:.95]); toc
%Elapsed time is 0.002386 seconds.

