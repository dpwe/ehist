function a = make_abbrev(s)
%  create an abbreviation for string or string array s
% 2014-01-03 Dan Ellis dpwe@ee.columbia.edu

% first letter plus any letters immediately following a space or
% underscore or .

if iscell(s)
  for i = 1:length(s)
    a{i} = make_abbrev(s{i});
  end
else
  s2 = [' ', s];
  spcs = (s2 == ' ' | s2 == '_' | s2 == '.');
  afterspc = find(spcs(1:end-1) & ~spcs(2:end));
  a = s(afterspc);
end
