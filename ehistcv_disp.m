function ehistcv_disp(CV, name, legend)
% ehistcv_disp(CV, name, legend)
%    Display a subband energy covariance in CV.
% 2013-11-05 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2;  name = ''; end
if nargin < 3;  legend = 0; end

if length(CV) == 0
  return
end

% assume mel
nfft = 512; sr = 8000; nfilts = size(CV,2); width = 2.0;
% just to get the frqs
[wts, frqs] = fft2melmx(nfft, sr, nfilts, width);

fix = 1:length(frqs);

imagesc(fix,fix,10*log10(max(1,CV))); axis xy
axis('square');

caxis([30 90])

title(name, 'interpreter', 'none');
%colorbar

ytl = str2num(get(gca, 'YTickLabel'));
set(gca, 'YTickLabel', round(frqs(ytl)));
ylabel('freq / Hz');

if legend
  %colorbar;
  xtl = str2num(get(gca, 'XTickLabel'));
  % Keep only first and last label
  xtl = xtl([1 end]);
  set(gca, 'XTick', xtl);
  set(gca, 'XTickLabel', round(frqs(xtl)));
  xlabel('freq / Hz');
else
  set(gca, 'XTick', []);
end
