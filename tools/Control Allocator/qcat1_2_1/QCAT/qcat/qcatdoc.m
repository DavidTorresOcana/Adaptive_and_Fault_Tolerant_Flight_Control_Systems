function qcatdoc

% QCATDOC - View QCAT documentation in Matlab's Help browser.
  
  % Get path
  s = which('qcatdoc');
  % Extract directory
  slash = strfind(s,'/');
  if isempty(slash)
    % PC platform -> look for backslash instead
    slash = strfind(s,'\');
  end
  s = s(1:slash(end-1));
  % Open documentation in browser.
  web(strcat(s,'doc/index.html'))
