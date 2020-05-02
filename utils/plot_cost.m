function plot_cost( name, plot_type )
%PLOT RESULTS: plot auxiliar function
% plot_results(opt_val, fn_values) produces some 
% plots of function outcome for each iteration
% and the difference errors used as a stopping criteria.


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

% check cost values were saved in a mat file
if exist(strcat(name,'/fx_k.mat'), 'file') ~= 2
    return;
end
load(strcat(name,'/fx_k.mat'))

if exist('x_true.mat', 'file') ~= 2 && ...
   exist('A.mat', 'file')      ~= 2 && ...
   exist('b', 'file') ~=2
    return;
end
load('x_true.mat'); load('A.mat'); load('b.mat');

if exist(strcat(name,'/opts.txt'), 'file') ~= 2
    return;
end
opts = table2struct(readtable(strcat(name,'/opts.txt')));
cost = lasso_function(A,x_true,b,opts.lambda);


% create directory to save png plots.
if ~exist('figures', 'dir')
    mkdir('figures');
end


% annotation
dim = [0.6, 0.7, 0.25, 0.15];
[step, lambda, maxiter, tol, ~ ] = parse_input_parameters( opts );
str = sprintf(['step = %g\n',...
               'lambda = %g\n',...
               'iterations = %d\n',...
               'tolerance = %.1g\n'],... 
               step,lambda,maxiter,tol);

% function to plot 1/k
gx = @(x) (fx_k(1)-cost).*x ;
k = 1./(1:length(fx_k)).^2;
l = 1./(1:length(fx_k));

x = 1000;
y = (fx_k(1)-cost) * k(x);


txt='$\frac{1}{k^2}$';
%------------------------------------------------------------
% plot minimization values with respect each iteration.
%------------------------------------------------------------
figname=sprintf('%s lambda=%.1e step=%.1e', upper(name),lambda, step);
figure('name', figname)
fn_vals = fx_k - cost;   
fn_plot{type}(fn_vals); 
hold on


h0=fn_plot{type}(gx(k));
h0.Color = [0.5 0.5 0.5]; h0.LineWidth = 0.5; h0.LineStyle = '--';
h1=fn_plot{type}(gx(l));
h1.Color = [0.5 0.5 0.5]; h1.LineWidth = 0.5; h1.LineStyle = '--';
text(x,y,txt,'interpreter','latex','HorizontalAlignment','right', 'FontSize',18)
hold off

%legend_names = strcat(upper(name));
%legend(legend_names,'FontSize',18,'Location','northeast');
annotation('textbox', dim, 'FontSize',14, 'String', str)

grid on
title(upper(name),'interpreter','none')
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontSize',18);
xlabel('${k}$','interpreter','latex','FontSize',18);

fig_name = sprintf('%s/%s_cost_%s.png',name,name,plot_type);
saveas(gcf,fig_name)

fig_name = sprintf('figures/%s_%s_lambda_%0.1g_step_%0.1g.png',name,plot_type,lambda, step);
saveas(gcf,fig_name);
end

