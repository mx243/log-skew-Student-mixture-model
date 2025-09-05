function m = sample_proGamma(a, b, D, typ)
% Usage: m = sample_proGamma(a, b, D, typ)
% This function uses 'ars.m' to draw a sample from proGamma(a, b, D, typ).

% Change Log:
%
%     1.1          05:jul:24    mx243      First version.
%     1.2          08:jul:24    mx243      Call error when isinf(a) && isinf (b), or a < - 30 * b, changed m when isinf(b / a).
%     1.14         23:aug:24    mx243      Added D. Now can sample from proGamma with D > 1.
%     1.16         03:sep:24    rfs34      Errors changed to warnings to avoid unexpected crashes.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.sample_proGamma.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

logf = @(x) log_density_proGamma(x, a, b, D, typ);
dlogf = @(x) deriv_log_density_proGamma(x, a, b, D, typ);

if(typ == 0)
    % Special cases first.
    if(isinf(a) && isinf(b))
        warning('sample_proGamma called with a = b = inf');
        m = 10; % Something to keep it going.
        return;
    elseif(isinf(a / (D * b)))
            if(a < 0)
                m = 1e100; % Probability is concentrated at a large value.
                return;
            end
            if(a > 0)
                m = 1e-100; % Probability is concentrated near 0.
                return;
            end
%    elseif(isinf(b / a))
%            m = 0.785; % Indeed it should be this value.
%            return;
    end 
    
    % Use ars.m to sample m. pts consists of points x where 
    % dlogf(x1) > 0, dlogf(x2) < 0.
    if(a + D * b <= 0)
        if(a / (D * b) < -30) % ARS performs poorly when a / (D * b) < -30.
            warning('sample_proGamma called with a < - 30 * D * b')
        end
        tmp = - a / (D * b);
        while(dlogf(exp(tmp)) <= 0)
            tmp = tmp - 1;
        end
        pts = [exp(tmp), exp(- a / (D * b)), 2 * exp(- a / (D * b))];
    else
        % Find lpos \in (0, 1] st. dlogf(lpos) > 0.
        const = 0;
        for j = 1 : D - 1
            const = const + log(1 + j / 2) - 1 / (2 + j);
        end
        const  = 2 * (a / (D * b) + 1 + const);
        if(const <= 0)
            lpos = 1;
        else
            lpos = min(1, 1 / const);
        end
        pts = [lpos, 3];
    end
    m = ars(logf, dlogf, pts, 0, inf);
end

if(typ == 1)
    % Special cases first.
    if(isinf(a) && isinf(b))
        warning('sample_proGamma called with a = b = inf');
        m = 10; % Something to keep it going.
        return;
    elseif(isinf(a / (D * b)))
            m = 1 + 1e-15; % Probability is concentrated near 1.
            return;
    elseif(isinf((D * b) / a))
            m = 1e100; % Probability is concentrated at a large value.
            return;
    end 
    
    % Use ars.m to sample m. pts consists of 2 points x1, x2, where 
    % dlogf(x1) > 0, dlogf(x2) < 0.
    const = - a / (D * b) - log(2 + (D - 1) / 2) + 1 / (4 + (D - 1));
    if(const >= 0)
        lpos = 2;
    else
        lpos = min(2, 1 + (exp(-1) - (D + 1) / 2) / const);
    end
    rpos = 1 + (D + 3) * D * b / (2 * a);
    pts = [lpos, rpos];
    m = ars(logf, dlogf, pts, 1, inf);
end

if(typ == 2)
    % Special cases first.
    if(isinf(a) && isinf(b))
        warning('sample_proGamma called with a = b = inf');
        m = 10; % Something to keep it going.
        return;
    elseif(isinf(a / b))
            m = 1e-100; % Probability is concentrated near 0.
            return;
    elseif(isinf(b / a))
            m = 1e100; % Probability is concentrated at a large value.
            return;
    end 

    % Use ars.m to sample m. pts consists of 2 points x1, x2, where 
    % dlogf(x1) > 0, dlogf(x2) < 0.
    pts = [b / (4 * a), (D + 1) * D * b / a];
    m = ars(logf, dlogf, pts, 0, inf);
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
