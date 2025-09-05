function y = log_total(x, dim)
% Usage: y = log_total(x, dim)
% Calculates log(sum(exp(x), dim)) with as little numerical
%  inaccuracy as easily feasible.

% Change Log:
%
%     1.1          28:sep:19    jc2062   As first received from JC2062.
%     1.2          10:oct:19    jc2062   As received from jc2062.
%     1.5          11:oct:19    rfs34    Dimension made explicit.
%     1.43         23:oct:19    rfs34    Comments only changed.

% Change Log (as part of mx243 code):
%
%     1.1          08:jul:24    mx243    Taken from above.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_total.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

m = max(x);
m(m == -Inf) = [0];
y2 = x - m;
y = log(sum(exp(y2), dim)) + m;

end

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
