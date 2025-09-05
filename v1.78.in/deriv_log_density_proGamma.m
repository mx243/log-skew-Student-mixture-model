function y = deriv_log_density_proGamma(m, a, b, D, typ)
% Usage: y = deriv_log_density_proGamma(m, a, b, D, typ)
% This function returns (d/dm)(log(f(m))) where m ~ proGamma(a, b, D, typ) 
% and f is proportional to P(m).

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.14         23:aug:24    mx243      Added D. 

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.deriv_log_density_proGamma.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if(typ == 0)
    y = - (a / b + D);
    for j = 0 : D - 1
        y = y - psi(m + j / 2);
    end
    y = y * b;
end

if(typ == 1)
    y = - (a / b + D) + D * log(m - 1) + D * (m + (D - 1) / 2) ./ (m - 1);
    for j = 0 : D - 1
        y = y - psi(m + j / 2);
    end
    y = y * b;
end

if(typ == 2)
    y = - (a / b + D) + D * log(m + (D - 1) / 2) + D;
    for j = 0 : D - 1
        y = y - psi(m + j / 2);
    end
    y = y * b;
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
