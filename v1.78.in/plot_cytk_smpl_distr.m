function plot_cytk_smpl_distr(cytk_prior, XLIM, Ndistr)
% Usage: plot_cytk_smpl_distr(cytk_prior, XLIM, Ndistr)
% This function generates Ndistr sample distributions of X in the logged space 
% from the cytokine priors. These should match the shape of demeaned,
% descaled log(cytokine) distributions.
% XLIM = [L, R] specifies the plotting range in the logged domain.

% Change Log:
%
%     1.1          18:mar:25    mx243      First version.
%     1.74         19:mar:25    rfs34      Version number harmonisation and fixed a spelling typo.
%     1.75         19:mar:25    mx243      Rewritten: now use the unchanged ASI priors to 
%                                          get sample distributions and then apply linear transformation.
%     1.76         21:mar:25    mx243      Rewritten: now only plot distributions of X.
%     1.77         22:mar:25    rfs34      Adjusted syntax to work with old Matlab also.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.plot_cytk_smpl_distr.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

L = XLIM(1); 
R = XLIM(2);
Npt = 1e3;
plotx = L : (R - L) / Npt : R;

[NaivePriorSamples, sd] = get_prior_smpl(cytk_prior, Ndistr);

LN0 = NaivePriorSamples.LN0;
LN1 = NaivePriorSamples.LN1;

Nskip = 1;
Phi = NaivePriorSamples.phi(LN0 : Nskip : LN1);
C = NaivePriorSamples.C(LN0 : Nskip : LN1);

% get_probs takes in plotting range in the non-logged space.
% log_Px is distribution in the logged domain, i.e. on plotx. It has size Ndistr * (Npt + 1). 
[log_Px, log_ddf, log_Px1bar] = get_probs(Phi, C, exp(plotx), 1, -1000); 

plot_all_curves(plotx, {exp(log_Px)}, 1, ones(1, size(log_Px, 1)), {'m-'}, 'prior samples for log(cytokine)', '', 'Pdf');

fig = gcf;
xlim(XLIM);

M = findobj(fig, 'Color', 'm', 'LineStyle', '-');
flag = 0;
while 1 
    figure(fig);
    if ~flag
        if isunix,
           fprintf('Press any key to continue or ctrl-c to exit\n');
           command = '';
        else
           command = input('Press enter to continue; Type in ''exit'' to terminate.', 's');
        end
    else
        if isunix,
           pause;
        else
           command = input('', 's');
        end
    end
    if ~strcmp(command, 'exit')
        visid = randi(Ndistr);
        for it_line = 1 : Ndistr
            if it_line == visid
                set(M(it_line), 'Visible', 'on');
            else
                set(M(it_line), 'Visible', 'off');
            end
        end
    else
        for it_line = 1 : Ndistr
            set(M(it_line), 'Visible', 'on');
        end
        break;
    end
    flag = 1;
end

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
