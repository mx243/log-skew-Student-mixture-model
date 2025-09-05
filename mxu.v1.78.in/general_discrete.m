function y = general_discrete(weights)
% Usage: y = general_discrete(weights)
% Draws a sample y from a finite discrete distribution a constant
%  multiple of whose probability values are given in weights.

% Change Log:
%
%     1.1          28:sep:19    jc2062   As first received from JC2062.
%     1.2          10:oct:19    jc2062   As received from jc2062.
%     1.21         16:oct:19    rfs34    Warning for all zero weights added.
%     1.28.1.4     22:oct:19    rfs34    Catches NaN and Inf weights.
%     1.39         22:oct:19    rfs34    Taken this code back into the main line.
%     1.43         23:oct:19    rfs34    Comments and precise warning circumstances 
%                                        and text only changed.
%     1.44         28:oct:19    rfs34    Trivial bug fixed.

% Change Log (as part of mx243 code):
%
%     1.1          08:jul:24    mx243    Taken from above, and modified by rfs34 to avoid using i as an iterator.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.general_discrete.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if any(isnan(weights)),
   warning('NaN weights will be treated as zero in general_discrete');
end

s = sum(weights(~isnan(weights)));

if isinf(s) && sum(isinf(weights)) > 1,
   warning('Multiple infinite weights being passed to general_discrete - will always return last such');
elseif s == 0,
   warning('All zero weights being passed to general_discrete - will always return 1');
end

v = s * rand;
j = 1;
t = weights(1);
while v > t
    j = j + 1;
    t = t + weights(j);
end
y = j;
end

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
