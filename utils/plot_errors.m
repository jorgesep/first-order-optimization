function plot_errors( name, plot_type )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%close all;

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

% check if error.mat file exists
if exist(strcat(name,'/errors.mat'), 'file') ~= 2
    return;
end

% load error matrix
load(strcat(name,'/errors.mat'))

[rows ~] = size(errors) ;


figure
for i=1:rows
    fn_plot{type}(errors(i,:));
    if i == 1 , hold on; end
end
hold off
grid on
title(name)

ylabel('error','FontSize',14);
xlabel('$x_{k}$','FontSize',14,'interpreter','latex');

legend({'$\frac{\|x_{k+1}-x_k\|_2}{\|x_k\|_2}$', ... 
        '$\frac{|x_{k+1}-x_k|_1}{\# x_k}$', ...
        '$\frac{|x_{k+1}-x_{true}|_1}{|x_{true}|_1}$', ...
        '$\frac{\|b-Ax\|_2}{\|b\|_2}$'}, ...
        'interpreter','latex', 'FontSize',18, 'Location', 'Best')

fig_name = sprintf('%s/errors_%s.png',name,plot_type);
saveas(gcf,fig_name)

end

