function plot_all( names, plot_type )
%PLOT RESULTS: plot auxiliar function
% plot_results(opt_val, fn_values) produces some 
% plots of function outcome for each iteration
% and the difference errors used as a stopping criteria.


% options for plot functions 
fn_plot = cell(4,1);
fn_plot(1) = {@(x)     plot(x,'LineWidth', 2.0) };
fn_plot(2) = {@(x) semilogx(x,'LineWidth', 2.0) };
fn_plot(3) = {@(x) semilogy(x,'LineWidth', 2.0) };
fn_plot(4) = {@(x)   loglog(x,'LineWidth', 2.0) };

% prepare plot type
keySet = {'plot','semilogx','semilogy','log'};
valueSet = [1 2 3 4];
M = containers.Map(keySet,valueSet);
type = M(plot_type);


% load matrices
if exist('x_true.mat', 'file') ~= 2 && ...
   exist('A.mat', 'file')      ~= 2 && ...
   exist('b', 'file') ~=2
    return;
end
load('x_true.mat'); load('A.mat'); load('b.mat');


% check cost values were saved in a mat file
fx = []; legend_names = {}; % empty array
for i=1:length(names)
    algorithm_name = names{i} ;
    if exist(strcat(algorithm_name,'/fx_k.mat'), 'file') == 2 && ...
       exist(strcat(algorithm_name,'/opts.txt'), 'file') == 2
   
        load(strcat(algorithm_name,'/fx_k.mat'))

        opts = table2struct(readtable(strcat(algorithm_name,'/opts.txt')));
        [h, lambda, ~, ~, ~ ] = parse_input_parameters( opts );
        fx = [fx; fx_k-lasso_function(A,x_true,b,lambda)]; %#ok<AGROW>

        legend_names{i} = sprintf('%s lambda = %g step=%g', upper(algorithm_name), lambda, h); %#ok<AGROW>
        clear fx_k; 
        clear opts;
    end
end


clear x_true; clear A; clear b;
if isempty(fx), return; end


% create directory to save png plots.
if ~exist('figures', 'dir')
    mkdir('figures');
end



% plot functions
[rows,~] = size(fx);
figure
for i=1:rows

    fn_plot{type}(fx(i,:));
    if i==1, hold on; end

end

% k = 1./(1:length(raw_vals)).^2;
% h=fn_plot{type}(raw_vals(1,1).*k);
% h.Color = uint8([17 17 17]); h.LineWidth = 0.5; h.LineStyle = '--';

hold off
grid on
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex','FontSize',16);
%set(get(gca,'ylabel'),'rotation',0)
xlabel('${k}$','interpreter','latex','FontSize',16);
title('Proximal-Gradient Algorithms','interpreter','none')
legend(legend_names,'FontSize',10,'Location','best', 'interpreter', 'none');


fig_name = sprintf('figures/functions_all_%s.png',plot_type);
saveas(gcf,fig_name);
end

