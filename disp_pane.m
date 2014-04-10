function disp_pane(X, CV, pane, maxpanes, name, tag)
% disp_pane(X, CV, pane, maxpanes, name, tag)
%   Fill in the information for a single pane
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; pane = 0; end
if nargin < 4; maxpanes = 0; end
if nargin < 5; name = ''; end
if nargin < 6; tag = ''; end

% Basic layout
rows = 4; 
cols = 2;
% maxpanes == rows*cols?

if maxpanes > 0 && pane > maxpanes
  % Don't try to plot off the end
  return
end

if pane > 0
  subplot(rows, 2 * cols, 2*pane-1);
end

ehist_disp(X, [], 0, name, (pane>(maxpanes-2)), mod(pane,2));
if length(tag)
  % add extra info on plot
  textypos = 0;
  textxpos = 5;
  text(textxpos, textypos, tag, ...
       'Color', [1 1 1], ...
       'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
end

% also cov
subplot(rows, 2*cols, 2*pane);
ehistcv_disp(CV, '', (pane>(maxpanes-2)));
