function [mytitl] = titlfunction(mytitl);
% Usage: [mytitl] = titlfunction(mytitl);
% This function manages titls, the record of which programs have been used
% in a given run.
% It expects the current titl of the caller to be in the parameter mytitl, 
% which should be persistent.
% It includes this in the global cell array titls, unless mytitl is empty; mytitl
% is returned as NaN so that it will not be included again.
%
% Typical usage may be found in the file ~rfs/matlab/typtitl .
%
% Global varible names used:
%   titls (for its standard purpose).


% Change Log:
%
%     1.1          25:may:01    rfs      First version saved in sccs.
%     1.2          25:may:01    rfs      Typical usage file referred to in comments.

% Change Log as part of mx243 code:
%
%     1.1          19:jun:24    mx243    First version saved in sccs.

%  /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.titlfunction.m 1.78 25/06/09 17:25:32


global titls
if ~isempty(mytitl) & ~isnan(mytitl),
   if isempty(titls),
      titls = cell(0, 1);
      titls{1} = ['Started ', datestr(now)];
   end
   fprintf('%s\n', mytitl);
   titls{length(titls) + 1} = mytitl;
   mytitl = NaN;
end




% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:


