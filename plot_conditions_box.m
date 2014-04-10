function plot_conditions_box(asr, titl, addlegend)
%  plot_conditions_box(asr, titl, addlegend)
%    plot box-and-whisker plot of ASR werr broken down by conditions
% 2014-01-06 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; titl = '<unk>'; end
if nargin < 3; addlegend = 0; end;

% plot box plot

gtot = max(asr.typeCode)+1;
nutts = length(asr.typeCode);

boxplot([asr.err./asr.chr, asr.err./asr.chr], ...
        [asr.typeCode, ...
         gtot*ones(1,nutts)]);

for i = 1:gtot-1
  N(i) = sum(asr.typeCode==i);
end
N(gtot) = sum(N);

set(gca, 'YGrid', 'on')
ymax = 1.0;
axis([0 gtot+1 0 ymax]);
set(gca, 'YTick', [0:0.2:1]);

set(gca, 'XTickLabelMode', 'auto');
set(gca, 'XTick', [1:gtot]);
set(gca, 'XTickLabel', [asr.typeNames, 'TOT']);

textypos = 0.02;
hold on; 
for i = 1:gtot
  text(i, textypos, num2str(N(i)), 'Color', [0 0 0], ...
       'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end 

hold off;

title([titl, ' - Error breakdown by condition'], 'interpreter', 'none');
xlabel('condition')
ylabel('error rate')

% Make whiskers easy to select in illustrator - unique color
h = [findobj(gca,'tag','Upper Whisker'), ...
     findobj(gca, 'tag', 'Lower Whisker')];
set(h, 'Color', [0.1 0.1 0.1]);
