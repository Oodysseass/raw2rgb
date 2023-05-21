close all
clc
set(0, 'DefaultFigureWindowStyle', 'docked');


% read .dng
[rawim, XYZ2Cam, wbcoeffs] = readdng("./RawImage.dng ");
[M0, N0] = size(rawim);


% % image in full scale and right bayer pattern 
bayertype = 'RGGB';
method = 'nearest';
[Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M0, N0);
figure(), imshow(Ccam);
figure(), imshow(Cxyz)
figure(), imshow(Clinear);
figure(), imshow(Csrgb);

% % histograms
figure()
subplot(2, 2, 1)
imhist(Ccam(:, :, 1));
subplot(2, 2, 2)
imhist(Ccam(:, :, 2));
subplot(2, 2, 3)
imhist(Ccam(:, :, 3));

figure()
subplot(2, 2, 1)
imhist(Cxyz(:, :, 1));
subplot(2, 2, 2)
imhist(Cxyz(:, :, 2));
subplot(2, 2, 3)
imhist(Cxyz(:, :, 3));

figure()
subplot(2, 2, 1)
imhist(Clinear(:, :, 1));
subplot(2, 2, 2)
imhist(Clinear(:, :, 2));
subplot(2, 2, 3)
imhist(Clinear(:, :, 3));

figure()
subplot(2, 2, 1)
imhist(Csrgb(:, :, 1));
subplot(2, 2, 2)
imhist(Csrgb(:, :, 1));
subplot(2, 2, 3)
imhist(Csrgb(:, :, 3));

% % demosaic using bilinear interpolation
method = 'linear';
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M0, N0);
figure(), imshow(Csrgb);

% % higher resolution
M = 7000;
N = 9000;
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

% % comparing images using smaller dimensions
M = 2500;
N = 3750;
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

M = 1238;
N = 1238;
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

% nearest with such small dimensions spoils the image
method = 'nearest';
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

% % wrong bayer patterns just for comparison
method = 'linear';
bayertype = 'BGGR';
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

bayertype = 'GBRG';
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);

bayertype = 'GRBG';
[Csrgb, ~, ~, ~] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
                                        bayertype, method, M, N);
figure(), imshow(Csrgb);