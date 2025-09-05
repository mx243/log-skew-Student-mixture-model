function test_calc_ASI(ASI_smpl_file)
% Usage: test_calc_ASI(ASI_smpl_file)
% This function plots the value of j(x1, x1_bar) for all k in a histogram, 
% then plots the curve of the distribution we obtained by training the Stu 
% model with the values of j(x1, x1_bar), and the mean and centiles of this
% distribution.

% Change Log:
%
%     1.1          18:sep:24    mx243      First version.
%     1.40         20:sep:24    rfs34      With edits and actions for mx243 after review by rfs34.
%     1.41         23:sep:24    mx243      Changes made following the instructions.
%     1.45         03:oct:24    rfs34      Merging changes made by rfs34 in the meantime;
%                                          now uses samples from LN0 to LN1, less any discard for convergence.
%     1.46         03:oct:24    rfs34      Removed fudged normalisation; changed way centiles are plotted.
%     1.50         09:oct:24    rfs34      And the way the mean is plotted, similarly.

mytitl = ' /home/rfs/ramakrishnan/software/lastact/mx243/SCCS/s.test_calc_ASI.m 1.78 25/06/09 17:25:32 ';

persistent mytitldone
if isempty(mytitldone),           
   mytitldone = titlfunction(mytitl);
end

global titls

load(ASI_smpl_file, 'Samples', 'x_ob', 'ASI', 'ASI025', 'ASI975', 'K', 'LN0', 'LN1');
Phi = Samples.phi(LN0 : LN1); 
C = Samples.C(LN0 : LN1);
smpl_num = LN1 - LN0 + 1;
p_dis = 0.15; % Discarding the first 15% of samples.

L = min(x_ob);
R = max(x_ob);
if L == R
    L = L - 4;
    R = R + 4;
end;
Nx1 = 100;
Nx2 = 1000;
x1 = L : (R - L) / Nx1 : R;
x2 = L : (R - L) / Nx2 : R;

log_y = -Inf(1, Nx2 + 1); 
cnt = 0;
for num = ceil(p_dis * smpl_num) : smpl_num  
    for it_c = 1 : C(num)
        tmp_log_y = log_density_skewStu(x2, Phi(num).mu(it_c), Phi(num).S(it_c), Phi(num).m(it_c), ...
                                        Phi(num).m(it_c) - 1, Phi(num).nu(it_c), 1);
        tmp_log_y = tmp_log_y + log(Phi(num).p(it_c));
        log_y = -Eadd(-log_y, -tmp_log_y);
    end
    cnt = cnt + 1;
    fprintf('\rDistr of j: ...done %d of %d...', cnt, smpl_num - ceil(p_dis * smpl_num) + 1);
end
log_y = log_y - log(cnt); 

% const = trapz(x2, exp(log_y)); 
const = 1; % It should already be normalised, and if it isn't then there's something wrong.
y = (exp(log_y) / const) * (K * (R - L) / Nx1); 

clf;
hold on;
hist(x_ob, x1);
plot(x2, y, 'g-', 'LineWidth', 1);
bins = hist(x_ob, x1);
H = max(bins);

% Plot the mean, with a line so that zooming in to it gets us a value rather than a barwidth.
plot([ASI, ASI], [0, H * 1.2], 'g--');

% Centiles are plotted as dashed lines.
plot([ASI025, ASI025], [0, H * 1.2], 'r--');
plot([ASI975, ASI975], [0, H * 1.2], 'r--');
hold off;

title('Distribution of ASI samples and mean posterior density');
xlabel('ASI in sample (nats)');

return;


% Local Variables: 
% indent-line-function: indent-relative
% eval: (auto-fill-mode 0)
% End:
