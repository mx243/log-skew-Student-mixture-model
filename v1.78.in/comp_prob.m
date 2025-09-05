function p = comp_prob(v1, v2)
% Usage: p = comp_prob(v1, v2)
% This function returns the mean posterior probability p of random 
% variable X being smaller than Y, where v1 is a vector of samples of X, 
% v2 is a vector of samples of Y, and the prior on p is the uniform
% distribution on [0, 1].

% Change Log:
%
%     1.1          17:mar:25    mx243      First version.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.comp_prob.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls


cnt = 0;
tot = size(v1, 1) * size(v2, 1);

for it_p1 = 1 : size(v1, 1)
    for it_p2 = 1 : size(v2, 1)
        if v1(it_p1) < v2(it_p2)
            cnt = cnt + 1;
        elseif v1(it_p1) == v2(it_p2)
            cnt = cnt + 1/2;
        end
    end
end

p = (cnt + 1) / (tot + 2);

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
