function plot_results(true_value, names, plot_type, varargin)
%PLOT RESULTS: plot auxiliar function
% plot_results(opt_val, fn_values) produces some 
% plots of function outcome for each iteration
% and the difference errors used as a stopping criteria.

close all;

% plot functions
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

% create directory to save png plots.
if ~exist('figures', 'dir')
    mkdir('figures');
end

% zeros will not be displayed on the graph.
%raw_vals = varargin{1}(i,:);
raw_vals = varargin{1};
raw_vals(raw_vals==0) = nan;

%------------------------------------------------------------
% plot minimization values with respect each iteration.
%------------------------------------------------------------

for i=1:length(names)
    
    fn_vals = raw_vals(i,1:end-1);
    fn_diff = fn_vals - raw_vals(i,2:end);
    
    fn_result  = fn_vals - true_value;
    
    figure
    fn_plot{type}(fn_vals);
    hold on;
    fn_plot{type}(fn_diff);
    fn_plot{type}(fn_result);
    hold off;
    
    ylabel('${f(x_k)}$','FontSize',14,'interpreter','latex');
    xlabel('$x_{k}$','FontSize',14,'interpreter','latex');
    
    str = sprintf('%s function convergence comparison alternatives', upper(names{i}));
    title(str);
    legend({'$f(x_k)$', ...
        '$f(x_{k-1})-f(x_k)$', ...
        '$f(x_{k})-f(x^{*})$'},...
        'FontSize',15,'interpreter','latex','Location','best');
    grid on
    
    fig_name = sprintf('figures/%s.png',names{i});
    saveas(gcf,fig_name)
end

%------------------------------------------------------------
% plot all functions in a single plot.
%------------------------------------------------------------
figure
for i=1:length(names)
    fn_vals = raw_vals(i,1:end) - true_value;   
    fn_plot{type}(fn_vals);
    if i == 1 , hold on; end
end
k = 1./(1:length(raw_vals)).^2;
h=fn_plot{type}(raw_vals(1,1).*k);
h.Color = uint8([17 17 17]); h.LineWidth = 0.5; h.LineStyle = '--';

hold off 
grid on
legend_names = strcat(upper(names'));
legend(legend_names,'FontSize',15,'Location','northeast');
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex');
xlabel('${k}$','interpreter','latex');
title('Proximal-Gradient Algorithms')


saveas(gcf,'figures/functions.png')
    

%------------------------------------------------------------
% plot errors.
%------------------------------------------------------------
for i=1:length(names)
    str = sprintf('%s_errors', names{i} );
    load( sprintf('%s.mat', str) );
    errors = eval(str);
    
    figure 
    fn_plot{type}(errors(1,:))
    hold on
    fn_plot{type}(errors(2,:))
    fn_plot{type}(errors(3,:))
    fn_plot{type}(errors(4,:))
    hold off
    grid on


    ylabel('error','FontSize',14);
    xlabel('$x_{k}$','FontSize',14,'interpreter','latex');
    
    str = sprintf('%s Error comparison alternatives', upper(names{i}));
    title(str);
    
    legend({'$\frac{\|x_{k+1}-x_k\|_2}{\|x_k\|_2}$', ... 
            '$\frac{|x_{k+1}-x_k|_1}{\# x_k}$', ...
            '$\frac{|x_{k+1}-x_{true}|_1}{|x_{true}|_1}$', ...
            '$\frac{\|b-Ax_k\|_2}{\|b\|_2}$'},...
            'interpreter','latex','FontSize',15,'Location','best')
    
    fig_name = sprintf('figures/%s_errors.png',names{i});
    saveas(gcf,fig_name)
    
    
end

%------------------------------------------------------------
% methods relative error comparison.
%------------------------------------------------------------
figure
for i=1:length(names)
    str = sprintf('%s_errors', names{i} );
    errors = eval(str);
    fn_plot{type}(errors(1,:));
    if i == 1 , hold on; end
    
end
hold off
grid on
ylabel('$\frac{\|x_{k+1}-x_k\|_2}{\|x_k\|_2}$','FontSize',16, ...
       'interpreter','latex')
set(get(gca,'ylabel'),'rotation',0)
xlabel('$x_k$','FontSize',14, 'interpreter','latex')
title('Norm Error convergence');
legend_names = strcat(upper(names'));
legend(legend_names,'FontSize',15,'Location','best');
saveas(gcf,'figures/norm_error.png')




end

