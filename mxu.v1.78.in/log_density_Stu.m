function y = log_density_Stu(x, mu, m, t, typ)
% Usage: y = log_density_Stu(x, mu, m, t, typ)
% This function returns log(P(x)) where x is a row vector of length K and 
% its entries are i.i.d ~ Stu(mu, S@1, m, r@r_from_t(m, t, typ)) 1D Student.
% This function is used when integrating wrt mu, m and t.

% Change Log:
%
%     1.1          31:jul:24    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.log_density_Stu.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

r = r_from_t(m, t, typ);
K = size(x, 2);

y = K * (- (1 / 2) * log(2 * pi) + gammaln(m + 1 / 2) - gammaln(m) + m .* log(r));
for j = 1 : K
    y = y - (m + 1 / 2) .* log(r + (1 / 2) * (x(j) - mu) .^ 2);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
