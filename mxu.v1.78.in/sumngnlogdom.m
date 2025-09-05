function [logout, logformula] = sumngnlogdom(N, logg);
% Usage: [logout, logformula] = sumngnlogdom(N, logg);
% Calculates the sum from j = 0 to infinity of j ^ N * g ^ j,
%  for abs(g) < 1, where g = exp(logg).
% N should be a non-negative integer and g a real scalar with abs(g) < 1.
% out is the numerical value for a particular g.
% logformula is a polynomial, with logformula(k) being the log of the coefficient of g^(k-1),
%  such that out = sum(exp(formula) .* g .^ [0 : N].') ./ (1 - g) .^ (N + 1).
% But this version of the function does everything in the log domain.

% Change Log:
%
%     1.1          24:aug:24    rfs34    First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sumngnlogdom.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if ~isscalar(N) || ~isreal(N) || N < 0 || N ~= round(N),
   error('N must be a non-negative integer');
end

if ~isreal(logg) || any(logg >= 0),
   error('g must be real and non-negative, and less than 1 for finite answer');
end

% We proceed by symbolic algebra,
% using the induction sumngn(N, g) = g * (d/dg)sumngn(N - 1, g),
% starting from sumngn(0, g) = 1 / (1 - g).

% We set the starting formula to 1, i.e. logformula to zero.

logformula = 0;

for n = 1 : N,

   % Take derivative of formula.
   logderivformula = log([1 : length(logformula) - 1].') + reshape(logformula(2 : end), [length(logformula) - 1, 1]);

   % Take n times formula.
   logdenom = log(n) + logformula;

   % Add the derivative.
   logdenom = - Eadd(- logdenom, - [logderivformula; -Inf]);

   % Subtract g times the derivative.
   logdenom = - Esub(- logdenom, - [-Inf; logderivformula]);

   % Multiply by g.
   logformula = [-Inf; logdenom];

end

logout = -Inf(size(logg));
on1 = ones(size(logg));
for n = N + 1 : -1 : 1,
   logformulan = logformula(n);
   logout = - Eadd(- (logout + logg), - logformulan(on1));
end
logout = logout - log(1 - exp(logg)) .* (N + 1);

ind = find(logg >= 0);
out(ind) = Inf;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
