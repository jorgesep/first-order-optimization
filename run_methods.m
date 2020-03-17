clear all;

% file names.
img_name = 'cameraman';
filename = strcat(img_name, '.jpg');
matfile = strcat('A_',img_name, '.mat');
eigfile = strcat('A_eigens_', img_name, '.mat');

% load image
img = im2double(imread(filename));

% vectorize true image
x = img(:);

% get blur kernel (gaussian) of window size nine and sigma 4.0
h = gaussian_kernel(9,4);

% linear operator
tic; fprintf('Getting linear operator ... ')
if exist(matfile, 'file') == 2
    S = load(matfile,'A');
    A = S.A;
else
    A = blur_optimize(h,img);
    save(matfile,'A');
end
elapsed=toc; fprintf('elapsed %f\n', elapsed)


% maximum eigen value
tic; fprintf('Computing eigenvalues ... ')
if exist(eigfile, 'file') == 2
    E = load(eigfile);
    A_square_eigs = E.A_square_eigs;
else
    A_square_eigs = eigs(A_square) ;
    save(eigfile,'A_square_eigs');
end
L = max(A_square_eigs);
elapsed = toc; fprintf('elapsed %f\n', elapsed)

% blurring the image
tic; fprintf('Blurring image ... ')
b = A * x;
elapsed=toc; fprintf('elapsed %f\n', elapsed)

% objective function
f_obj = @(A, x, b, lambda) norm(b-A*x)^2 + lambda*sum(abs(x));

% parameters for ISTA 
x_0 = b;  % init with blurred image
step = 1 / L;
lambda = 2e-5;
maxiter = 1000;

step = step * 0.5;
[x_ista, f_ista] = ista(A, b, x_0, step, lambda, maxiter, f_obj);

[x_fista, f_fista] = fista(A, b, x_0, step, lambda, maxiter, f_obj);

%step = lambda;
%[x_grad, f_grad] = gradient_descent(A, b, x_0, step, lambda, maxiter, f_obj);

fx_min = f_obj(A,x,b,lambda);
figure(1);
semilogx(f_ista); hold on;
semilogx(f_fista); 
%semilogx(f_grad);
hold off;

ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('Algoritmos Gradiente-Proximal','interpreter','latex', 'FontWeight','bold')
legend({'ISTA','FISTA'},'Location','northeast')
grid on;

figure(2); 
loglog(f_ista); hold on;
loglog(f_fista); 
%loglog(f_grad);
hold off;
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('Algoritmos Gradiente-Proximal','interpreter','latex', 'FontWeight','bold')
legend({'ISTA','FISTA'},'Location','northeast')
grid on;

figure(3);
vals_ista = f_ista-fx_min ;
vals_fista = f_fista-fx_min ;
%vals_grad = f_grad-fx_min ;
semilogx(vals_ista(vals_ista>0)); hold on;
semilogx(vals_fista(vals_fista>0)); 
%semilogx(vals_grad(vals_grad>0));
hold off;
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('Algoritmos Gradiente-Proximal','interpreter','latex', 'FontWeight','bold')
legend({'ISTA','FISTA'},'Location','northeast')
grid on;

figure(4);
loglog(vals_ista(vals_ista>0)); hold on;
loglog(vals_fista(vals_fista>0)); 
%loglog(vals_grad(vals_grad>0)); 
hold off;
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('Algoritmos Gradiente-Proximal','interpreter','latex', 'FontWeight','bold')
legend({'ISTA','FISTA'},'Location','northeast')
grid on;

figure(5);
semilogx(f_ista(1:end-1)-f_ista(2:end)); hold on;
semilogx(f_fista(1:end-1)-f_fista(2:end)); 
%semilogx(f_grad(1:end-1)-f_grad(2:end));
hold off;
ylabel('${f(x_k)-f(x^{*})}$','interpreter','latex', 'FontWeight','bold')
xlabel('${k}$','interpreter','latex', 'FontWeight','bold')
title('Algoritmos Gradiente-Proximal','interpreter','latex', 'FontWeight','bold')
legend({'ISTA','FISTA'},'Location','northeast')
grid on;

% % reshape the blurred image
img_blurred = reshape(b,size(img)); 
img_ista    = reshape(x_ista,size(img));
img_fista   = reshape(x_fista,size(img));
% 
% % Display images
% figure(6); imshow(im2uint8(img_blurred)); title('Blurred');
% figure(7); imshow(im2uint8(img_ista));title('ISTA');
% figure(8); imshow(im2uint8(img_fista));title('FISTA');
% 
% figure(4)
% plot(fopt(fopt>0));
iptsetpref('ImshowBorder','tight');
figure(6);
subplot(1,3,1), imshow(im2uint8(img_blurred)); title('Blurred');
subplot(1,3,2), imshow(im2uint8(img_ista));title('ISTA');
subplot(1,3,3), imshow(im2uint8(img_fista));title('FISTA');
%tightfig;
%spaceplots([0 0 0 0], [.02 .02]);
%set(gca, 'LooseInset', get(gca,'TightInset'))
% F= getframe(f);
% img1 = F.cdata;
% imwrite(img1, 'test.jpg');

% %Matlab runs a script called ?startup.m? when it starts. This is located in your ~/matlab/ 
% %remove the gray border
% iptsetpref('ImshowBorder','tight');
% %removes menu and toolbar from all new figures
% set(0,'DefaultFigureMenu','none');
% %makes disp() calls show things without empty lines
% format compact;
% set(0,'Default');
Fvalues=[f_ista-fx_min; f_fista-fx_min];
Xvalues=[x_ista; x_fista];
save('Fvalues.mat','Fvalues');
save('Xvalues.mat','Xvalues');
