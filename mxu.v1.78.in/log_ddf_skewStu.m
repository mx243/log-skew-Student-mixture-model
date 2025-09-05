function v = log_ddf_skewStu(x_ob, mu, S, m, r, nu, D, nalpha)
% Usage: v = log_ddf_skewStu(x_ob, mu, S, m, r, nu, D, nalpha)
% This function returns the log of the decumulative probability 
% P(x1 > x_ob(1), x1_bar @ x_ob(2 : D)), where x ~ skewStu(mu, S, m, r, nu, D). 
% The input x can be a D * K matrix and the output will be a 1 * K row vector. 
% If nalpha > 0, then the probability P(x1 > d) is computed using Monte-Carlo integration over
%  alpha, the number of samples of alpha drawn is specified by nalpha.
% If nalpha < 0, then numerical integration over alpha is used using abs(nalpha) + 1 points.
%

% Change Log:
%
%     1.1          25:sep:24    mx243      First version.
%     1.45         03:oct:24    rfs34      Layout only changed.
%     1.61         26:oct:24    rfs34      Bug in use of ones fixed.
%     1.66         01:nov:24    mx243      Minor changes to increase speed.
%     1.67         03:nov:24    mx243      Increased vectorisation to increase speed. 
%     1.71         25:jan:25    rfs34      Added option to numerically integrate over log-uniformly 
%                                          spaced values of alpha, but not yet debugged.
%     1.72         25:jan:25    rfs34      Fixed sign error in above.
%     1.75         20:mar:25    mx243      Fixed trivial bug on line 38 (in v1.74).

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_ddf_skewStu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

oldmatlab = verLessThan('matlab', '8.0.1');

K = size(x_ob, 2);

usenumint = 0;
if nalpha > 0,
   alpha = sample_Gamma(m, r, [nalpha, 1]);
   EGamma = 0;
elseif nalpha < 0,
   usenumint = 1;
   nalpha = abs(nalpha) + 1;
   logalphalow = log(1e-4 * m / r);
   logalphahigh = log(1e2 * m / r);
   deltaalpha = (logalphahigh - logalphalow) / (nalpha - 1);
   alpha = exp([logalphalow : deltaalpha : logalphahigh].');
   EGamma = - log_density_Gamma(alpha, m, r) - log(alpha);
else
   error('nalpha may not be zero in log_ddf_skewStu.m');
end

sqrtalpha = sqrt(alpha);

if oldmatlab,
   y = x_ob - mu(:, ones(1, K));
else
   y = x_ob - mu; % Same as x_ob - mu(:, ones(1, K)) but faster on my machine.
end
y1bar = y(2 : D, :);
nu1bar = nu(2 : D, 1);
S_bar = S(2 : D, 2 : D);
S1 = S(1, 2 : D);

Sy = S_bar * y1bar;
ySy = sum(y1bar .* Sy, 1);
nuSy = nu1bar' * Sy;
S1y = S1 * y1bar;

S1nu = S1 * nu1bar;
nuSnu = nu1bar' * (S_bar * nu1bar);

tmpx = sqrt(S(1, 1) / 2) * (sqrtalpha * (y(1, :) + S1y / S(1, 1)) - nu(1) - S1nu / S(1, 1));
const = (1 / 2) * (log(pi) - log(2) - log(S(1, 1)) - D * log(2 * pi) - Edet(S) + S1nu ^ 2 / S(1, 1) - nuSnu);

% v has size [nalpha, K] at this point.
v = (alpha / 2) * (S1y .^ 2 / S(1, 1) - ySy) + sqrtalpha * (nuSy - S1y * (S1nu / S(1, 1)));

term = ((D - 1) / 2) * log(alpha) - EGamma;
if oldmatlab,
   v = v + term(:, ones(1, K));
else
   v = v + term;
end

v = v + log(erfc(tmpx));

% v has size [1, K] now.
v = -Esum(-v, 1) + const;
if usenumint,
   v = v + log(deltaalpha);
else
   v = v - log(nalpha);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
