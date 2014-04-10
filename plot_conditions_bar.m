function plot_conditions_bar(asr, titl, addlegend)
%  plot_conditions_bar(asr, titl, addlegend)
%    plot bar chart broken down by conditions
% 2014-01-03 Dan Ellis dpwe@ee.columbia.edu

if nargin < 2; titl = '<unk>'; end
if nargin < 3; addlegend = 0; end;

% Copy across the names of each type
cond.names = asr.typeNames;

% sum up by conditions
ncond = length(cond.names);
for i = 1:ncond
  ix = find(asr.typeCode == i);
  %disp(dem.types.names(i));
  %for ii = 1:min(5,length(ix)); disp(asr.uttID{ix(ii)}); end

  cond.chr(i) = sum(asr.chr(ix));
  cond.snt(i) = sum(asr.snt(ix));
  cond.corr(i) = sum(asr.corr(ix));
  cond.sub(i) = sum(asr.sub(ix));
  cond.del(i) = sum(asr.del(ix));
  cond.ins(i) = sum(asr.ins(ix));
  cond.err(i) = sum(asr.err(ix));
  cond.serr(i) = sum(asr.serr(ix));
  cond.N(i) = length(ix);
end

% where to put grand totals
gtot = ncond + 1;

cond.chr(gtot) = sum(cond.chr);
cond.snt(gtot) = sum(cond.snt);
cond.corr(gtot) = sum(cond.corr);
cond.sub(gtot) = sum(cond.sub);
cond.del(gtot) = sum(cond.del);
cond.ins(gtot) = sum(cond.ins);
cond.err(gtot) = sum(cond.err);
cond.serr(gtot) = sum(cond.serr);
cond.N(gtot) = sum(cond.N);

% plot bar chart

bar(1:gtot, ([cond.sub; cond.del; cond.ins]*diag(1./cond.chr))', 0.5, 'stack');

set(gca, 'YGrid', 'on')
ymax = 0.7;
axis([0 gtot+1 0 ymax]);

for i = 1:ncond
%  blabel{i} = [cond.names{i}, '', num2str(cond.N(i))];
  blabel{i} = cond.names{i};
end
%blabel{gtot} = ['TOT',num2str(cond.N(gtot))];
blabel{gtot} = 'TOT';

textypos = 0.0;
hold on; 
for i = 1:gtot
  text(i, textypos, num2str(cond.N(i)), 'Color', [1 1 1], ...
       'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
end 
hold off;

set(gca, 'XTickLabel', blabel);
title([titl, ' - Error breakdown by condition'], 'interpreter', 'none');
xlabel('condition')
ylabel('error rate')

if addlegend
  legend('sub','del','ins',4)
end

