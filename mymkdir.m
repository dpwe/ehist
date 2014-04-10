function r = mymkdir(dir)
% r = mymkdir(dir)
%   Ensure that dir exists by creating all its parents as needed.
% 2006-08-06 dpwe@ee.columbia.edu

r = 0;

if length(dir) == 0
  return;
end

[x,m,i] = fileattrib(dir);
if x == 0
  [pdir,nn] = stripdir(dir);
  disp(['creating ',dir,' ... ']);
  r = mymkdir(pdir);
  % trailing slash results in empty nn
  if length(nn) > 0
    if length(pdir) == 0
      pdir = pwd;
    end
    mkdir(pdir, nn);
    r = 1;
  end
end

%%%%%%%%%%%%%%%%
function [head,tail] = stripdir(dirname)
% Divide off the last tail part of a dir path
lastslash = max([0,find(dirname == '/')]);
tail = dirname(lastslash+1:end);
head = dirname(1:lastslash-1);  % will return nothing when lastslash == 0
