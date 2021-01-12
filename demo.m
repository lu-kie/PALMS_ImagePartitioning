clear;
addpath(genpath('src'))
% Get example image (picture credits: Tuxyso / Wikimedia Commons)
f = im2double(imread('redMacaw.jpg'));

%%
% set jump penalty (larger choice -> less segments)
gamma = 0.75;
% Perform affine-linear image partitioning
[u,partition,a,b,c] = affineLinearPartitioning(f,'gamma',gamma);
% u: piecewise affine-linear approximation of f
% partition: corresponding partition
% a,b,c: slopes and offset (in matrix origin) of the first order
% polynomials in each pixel

% Show results
boundaries = boundarymask(partition);

figure,
subplot(2,2,1),imshow(f),title('Input image')
subplot(2,2,2),imshow(u),title('Pcw. affine approximation')
subplot(2,2,3),imshow(label2rgb(partition,'jet')),title('Partitioning')
subplot(2,2,4),imshow(1-boundaries),title('Boundaries')
