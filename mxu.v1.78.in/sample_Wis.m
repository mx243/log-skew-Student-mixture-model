function S = sample_Wis(m, R, D)
% Usage: S = sample_Wis(m, R, D)
% This function returns a sample from S ~ Wis(m, R, D).

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.2          15:jul:24    mx243      Draw multiple samples from Gamma and Normal
%                                          at once to increase speed.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_Wis.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

if(0)
C = zeros(D, D);
for row = 1 : D
    for col = row : D
        if(row == col)
            C(row, row) = sqrt(sample_Gamma(m + (D - row) / 2, 1));
        else
            C(row, col) = (1 / sqrt(2)) * randn;
        end
    end
end
end

tmp_m_vec = m + (D - 1) / 2 : - (1 / 2) : m;
tmp_Gam = sample_Gamma(tmp_m_vec, 1);
C = diag(sqrt(tmp_Gam));

tmp_Norm = randn([1, (D ^ 2 - D) / 2]);
cnt = 1;
for row = 1 : D - 1
    C(row, row + 1 : D) = (1 / sqrt(2)) * tmp_Norm(cnt : cnt + D - row - 1);
    cnt = cnt + D - row;
end

W = forcechol(inv(R));
C = C * W;
S = C' * C;

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
