function PP = ehist_disp(X, pcntls, maxf, name, xlegend, ylegend)
% PP = ehist_disp(X, pcntls, maxf, xlegend, ylegend)
%    Display an energy histogram encoded in X.
%    Returns actual percentile thresholds.
% 2013-11-05 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  pcntls = []; end
if nargin < 3;  maxf = 0; end
if nargin < 4;  name = ''; end
if nargin < 5;  xlegend = 0; end
if nargin < 6;  ylegend = 0; end

if length(pcntls) == 0 
  pcntls = [.05 .95];
end

if length(X) == 0
  PP = [];
  return
end

% Normalize by total count (assume same for each column)
sx = sum(X(:,1));
X = X/sx;

PP = histpercentile(X, pcntls);

if maxf > 0
  hzperbin = maxf/size(X, 2);
  frqs = [0:size(X, 2)-1] * hzperbin;
else
  % assume mel
  nfft = 512; sr = 8000; nfilts = size(X,2); width = 2.0;
  % just to get the frqs
  [wts, frqs] = fft2melmx(nfft, sr, nfilts, width);
end

fix = 1:length(frqs);
dbperbin = 1.0;

dd = [0:size(X, 1)-1] * dbperbin;

imagesc(fix,dd,X); axis xy
% Common axis
caxis([0 0.05]);

if ylegend
  ylabel('level / dB');
end

title(name, 'interpreter', 'none');
%colorbar

if size(PP,1) > 0
  hold on
  for i = 1:size(PP,1)
    %if (i==1); c = 'k'; else c = 'w'; end
    c = 'w';
    plot(fix, PP(i,:),c);
  end
  hold off
end

if xlegend
  xtl = str2num(get(gca, 'XTickLabel'));
  set(gca, 'XTickLabel', round(frqs(xtl)));
  xlabel('freq / Hz');
else
  set(gca, 'XTick', []);
end

