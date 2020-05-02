function  plot_x1x2( name )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



% read options from txt file
opts = inf;
opts_file = strcat(name,'/opts.txt') ;
if exist(opts_file, 'file') == 2
    opts = table2struct(readtable(opts_file));
end

ix1 = inf; ix2 = inf ;
if isfield(opts, 'ix1'), ix1 = opts.ix1; end
if isfield(opts, 'ix2'), ix2 = opts.ix2; end

if exist('x_true.mat', 'file') == 2
    load('x_true.mat');
end

% load x1 x2 array
if exist(strcat(name,'/x1x2.mat'), 'file') ~= 2
    return;
end

load(strcat(name,'/x1x2.mat'))


figure
plot(x1x2(1,:),x1x2(2,:),'LineWidth', 2.0); %#ok<NODEF>
hold on
plot(x1x2(1,1),x1x2(2,1),'g*','LineWidth', 2.0);
if ix1 ~= inf
    plot(x_true(ix1),x_true(ix2), 'rx', 'LineWidth', 2.0)
end
plot([x1x2(1,1) x_true(ix1)],[x1x2(2,1) x_true(ix2)],'--','Color', [0.5 0.5 0.5], 'LineWidth', 0.5);

hold off
grid on
tname = sprintf('%s convergence of two variables', upper(name));
title(tname,'interpreter','none','FontSize',16)

ylabel('$x_2$','FontSize',18,'interpreter','latex');
xlabel('$x_1$','FontSize',18,'interpreter','latex');

% legend({'$\frac{\|x_{k+1}-x_k\|_2}{\|x_k\|_2}$', ... 
%         '$\frac{|x_{k+1}-x_k|_1}{\# x_k}$', ...
%         '$\frac{|x_{k+1}-x_{true}|_1}{|x_{true}|_1}$', ...
%         '$\frac{\|b-Ax\|_2}{\|b\|_2}$'}, ...
%         'interpreter','latex', 'FontSize',18, 'Location', 'Best')

% anotation box
dim = [0.7, 0.3, 0.0, 0.0];
str = sprintf('vector position:\n$x_1$=%d\n$x_2$=%d\n', ix1, ix2);
%annotation('textbox', dim, 'FontSize',14, 'String', str,'interpreter','latex', 'FitBoxToText','on')           
annotation('textbox', dim, 'FontSize',14, 'String', str,'interpreter','latex')           
 
% ylim=get(gca,'ylim');
% xlim=get(gca,'xlim');
% text(xlim(1),ylim(2),str, 'FontSize',14,'interpreter','latex');

fig_name = sprintf('%s/x1x2_%s.png',name);
saveas(gcf,fig_name)

end

