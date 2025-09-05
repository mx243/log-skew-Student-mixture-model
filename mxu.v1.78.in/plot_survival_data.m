function plot_survival_data(data, censored, xr, colstring);
% Usage: plot_survival_data(data, censored, xr, colstring);
% Plots Kaplan-Meier-style plots of the deaths that occurred.
% data should be a vector of positive reals giving times of
%  either death or censoring. 
% censored should be a logical array of the same size, 
%  1 for patients for whom data is a censoring time,
%  0 for patients for whom data is a death time.

% Change Log:
%
%     1.1          28:sep:19    jc2062   As first received from JC2062.
%     1.2          10:oct:19    jc2062   As received from jc2062.
%     1.7          12:oct:19    rfs34    Changed parameters.
%     1.14         14:oct:19    rfs34    Accepted censored parameter.
%     1.34         19:oct:19    rfs34    Removed downward bias from calculation.
%     1.43         23:oct:19    rfs34    Comments only changed.
%     1.48         09:nov:19    rfs34    Switched ' to (:) to compensate for bug in hist: 
%                                        nature of output vector when input data empty.
%     1.49         12:nov:19    rfs34    Also need to cope with tt having length 1 for hist,
%                                        and removed unused return parameter.

% Change Log as part of mx243 code:
%
%     1.1          11:oct:24    mx243    Set Linewidth = 1;

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_survival_data.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

    N = length(data);
    [data, ind] = sort(data);
    censored = censored(ind);

    tt = unique(data);
    tt = tt(~isinf(tt));

    if length(tt) > 1,
       actuallydead = data(~censored & ~isinf(data));
       Ndeadattt = hist(actuallydead, tt);
       Ndeadattt = Ndeadattt(:); % hist doesn't reliably return the right sort of vector
                                 % when data is empty.
       Ndeadorcensoredattt = hist(data(~isinf(data)), tt);
       Ndeadorcensoredattt = Ndeadorcensoredattt(:);
    elseif length(tt) == 1,
       Ndeadattt = sum(~censored);
       Ndeadorcensoredattt = length(censored);
    else % i.e. if tt is empty
       Ndeadattt = zeros(0, 1);
       Ndeadorcensoredattt = zeros(0, 1);
    end

    Naliveatttminus = N - [0; cumsum(Ndeadorcensoredattt(1 : end - 1))];

    fracdieattt = 1 - Ndeadattt ./ Naliveatttminus;

    deadfrac = cumprod(fracdieattt);

    plottt = [tt.'; tt.'];
    plottt = plottt(:);
    plottt = [0; plottt];

    plotdeadfrac = [deadfrac.'; deadfrac.'];
    plotdeadfrac = plotdeadfrac(:);
    plotdeadfrac = [1; 1; plotdeadfrac(1 : end - 1)];

    plot(plottt, plotdeadfrac, colstring, 'LineWidth', 1);
    axis([0, xr, -0.01, 1.01]);

end

% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
