function [Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, ...
                                        wbcoeffs, bayertype, method, M, N)
    rawim = resize(rawim, M, N);
    
    if bayertype == "RGGB"
        Ccam = rggb(rawim, wbcoeffs, method, M, N);
    elseif bayertype == "BGGR"
        Ccam = bggr(rawim, wbcoeffs, method, M, N);
    elseif bayertype == "GBRG"
        Ccam = gbrg(rawim, wbcoeffs, method, M, N);
    elseif bayertype == "GRBG"
        Ccam = grbg(rawim, wbcoeffs, method, M, N);
    else
        Csrgb = [];
        Clinear = [];
        Cxyz = [];
        Ccam = [];
        fprintf("Not compatible bayer pattern\n")
        return
    end

    % % initialize transformation arrays
    cam2xyz = XYZ2Cam ^ -1;
    xyz2rgb = [3.2406 -1.5372 -0.4986;
                -0.9689 1.8758 0.0415;
                0.0557 -0.2040 1.0570];
    % normalize rows' sum to 1
    cam2xyz = cam2xyz ./ repmat(sum(cam2xyz, 2), 1, 3);
    xyz2rgb = xyz2rgb ./ repmat(sum(xyz2rgb, 2), 1, 3);

    % calculate each layer of Cxyz
    r = cam2xyz(1, 1) * Ccam(:,:,1) + cam2xyz(1, 2) * Ccam(:,:,2) + ...
        cam2xyz(1, 3) * Ccam(:,:,3);
    g = cam2xyz(2, 1) * Ccam(:,:,1) + cam2xyz(2, 2) * Ccam(:,:,2) + ...
        cam2xyz(2, 3) * Ccam(:,:,3);
    b = cam2xyz(3, 1) * Ccam(:,:,1) + cam2xyz(3, 2) * Ccam(:,:,2) + ...
        cam2xyz(3, 3) * Ccam(:,:,3);
    % concatenate arrays
    Cxyz = cat(3, r, g, b);
    Cxyz = max(0, min(1, Cxyz));

    % calculate each layer of Clinear
    r = xyz2rgb(1, 1) * Cxyz(:,:,1) + xyz2rgb(1, 2) * Cxyz(:,:,2) + ...
        xyz2rgb(1, 3) * Cxyz(:,:,3);
    g = xyz2rgb(2, 1) * Cxyz(:,:,1) + xyz2rgb(2, 2) * Cxyz(:,:,2) + ...
        xyz2rgb(2, 3) * Cxyz(:,:,3);
    b = xyz2rgb(3, 1) * Cxyz(:,:,1) + xyz2rgb(3, 2) * Cxyz(:,:,2) + ...
        xyz2rgb(3, 3) * Cxyz(:,:,3);
    % concatenate arrays
    Clinear = cat(3, r, g, b);
    Clinear = max(0, min(1, Clinear));
    % fix brightness
    grayim = rgb2gray(Clinear);
    grayscale = 0.25 / mean(grayim(:));
    Clinear = min(1, Clinear * grayscale);

    % Csrgb with approximation
    Csrgb = Clinear .^ (1 / 2.2);
end


function sampled = resize(rawim, M, N)
    % M0xN0 size of image
    [M0, N0] = size(rawim);

    % step of sampling rate in each dimension
    stepM = M0 / M * 2;
    stepN = N0 / N * 2;

    sampled = zeros(M, N);

    % % sample each color
    % two ways because when M0 / M and N0 / N is integer it is faster
    if mod(M0, M) == 0 && mod(N0, N) == 0
        % reverse order of sampling to make sure that
        % the edges of the grid is literally the same value
        sampled(1:2:end, 1:2:end) = rawim(1:stepM:end, 1:stepN:end);
        sampled(1:2:end, end:-2:2) = rawim(1:stepM:end, end:-stepN:2);
        sampled(end:-2:2, 1:2:end) = rawim(end:-stepM:2, 1:stepN:end);
        sampled(2:2:end, end:-2:2) = rawim(2:stepM:end, end:-stepN:2);
    else
    % second way is the classic with 2 fors
    % it essentially does the same thing but with an decimal step
    % which is rounded to the nearest censor of the same color
    % ii and jj are essentially dummy variables, but they help
    % for clarity
        % first censors that are on odd rows and cols
        ii = 1;
        tempi = 1;
        % take edge value first
        sampled(1, 1) = rawim(1, 1);
        for i = 1 : 2 : M
            jj = 1;
            for j = 3 : 2 : N
                jj = jj + stepN;
                % different rounding tactics so that
                % we don't need edge cases
                if N < N0
                    tempj = 2 * round(jj / 2) + 1;
                else
                    tempj = 2 * round(jj / 2) - 1;
                end
                sampled(i, j) = rawim(tempi, tempj);
            end
            ii = ii + stepM;
            if M < M0
                tempi = 2 * round(ii / 2) + 1;
            else
                tempi = 2 * round(ii / 2) - 1;
            end
        end
    
        % censors on odd rows and even cols
        ii = 1;
        tempi = 1;
        sampled(1, end) = rawim(1, end);
        for i = 1 : 2 : M
            jj = N0;
            for j = N - 2 : -2 : 2
                jj = jj - stepN;
                tempj = round(jj);
                if mod(tempj, 2) == 1
                    tempj = tempj + 1;
                end
                sampled(i, j) = rawim(tempi, tempj);
            end
            ii = ii + stepM;
            if M < M0
                tempi = 2 * round(ii / 2) + 1;
            else
                tempi = 2 * round(ii / 2) - 1;
            end
        end
    
        % censors on even rows and odd cols
        ii = M0;
        tempi = M0;
        sampled(end, 1) = rawim(end, 1);
        for i = M : -2 : 2
            jj = 1;
            for j = 3 : 2 : N
                jj = jj + stepN;
                if N < N0
                    tempj = 2 * round(jj / 2) + 1;
                else
                    tempj = 2 * round(jj / 2) - 1;
                end
                sampled(i, j) = rawim(tempi, tempj);
            end
            ii = ii - stepM;
            tempi = round(ii);
            if mod(tempi, 2) == 1
                tempi = tempi + 1;
            end
        end
    
        % censors on even rows and even cols
        ii = M0;
        tempi = M0;
        sampled(end, end) = rawim(end, end);
        for i = M : -2 : 2
            jj = N0;
            for j = N : -2 : 2
                jj = jj - stepN;
                tempj = round(jj);
                if mod(tempj, 2) == 1
                    tempj = tempj + 1;
                end
                if tempj <= 0
                    tempj = 2;
                end
                sampled(i, j) = rawim(tempi, tempj);
            end
            ii = ii - stepM;
            tempi = round(ii);
            if mod(tempi, 2) == 1
                tempi = tempi + 1;
            end
        end
    end
    
end


function Ccam = rggb(rawim, wbcoeffs, method, M, N)
    % % extract color components
    red = zeros(M, N);
    blue = zeros(M, N);
    red(1:2:end, 1:2:end) = rawim(1:2:end, 1:2:end);
    blue(2:2:end, 2:2:end) = rawim(2:2:end, 2:2:end);
    green = rawim - red - blue;

    % % white balance
    % mask with the size of image
    % leave green as it is
    wbmask = zeros(M, N);
    wbmask(1:2:end, 1:2:end) = wbcoeffs(1);
    wbmask(2:2:end, 2:2:end) = wbcoeffs(3);
    red = red .* wbmask;
    blue = blue .* wbmask;

    % interpolation
    if method == "nearest"
        Ccam = nearestRGGB(red, green, blue);
    elseif method == "linear"
        Ccam = bilinear(red, green, blue);
    else
        fprintf("Method not supported\n")
        Ccam =[];
    end
end


function Ccam = bggr(rawim, wbcoeffs, method, M, N)
    red = zeros(M, N);
    blue = zeros(M, N);
    red(2:2:end, 2:2:end) = rawim(2:2:end, 2:2:end);
    blue(1:2:end, 1:2:end) = rawim(1:2:end, 1:2:end);
    green = rawim - red - blue;

    wbmask = zeros(M, N);
    wbmask(2:2:end, 2:2:end) = wbcoeffs(1);
    wbmask(1:2:end, 1:2:end) = wbcoeffs(3);
    red = red .* wbmask;
    blue = blue .* wbmask;

    if method == "nearest"
        % bggr is essentially the same as rggb with r and b flipped
        Ccam = nearestRGGB(blue, green, red);

        red = Ccam(:, :, 3);
        Ccam(:, :, 3) = Ccam(:, :, 1);
        Ccam(:, :, 1) = red;
    elseif method == "linear"
        Ccam = bilinear(red, green, blue);
    else
        fprintf("Method not supported\n")
        Ccam =[];
    end
end


function Ccam = gbrg(rawim, wbcoeffs, method, M, N)
    red = zeros(M, N);
    blue = zeros(M, N);
    red(2:2:end, 1:2:end) = rawim(2:2:end, 1:2:end);
    blue(1:2:end, 2:2:end) = rawim(1:2:end, 2:2:end);
    green = rawim - red - blue;

    wbmask = zeros(M, N);
    wbmask(2:2:end, 1:2:end) = wbcoeffs(1);
    wbmask(1:2:end, 2:2:end) = wbcoeffs(3);
    red = red .* wbmask;
    blue = blue .* wbmask;

    if method == "nearest"
        Ccam = nearestGBRG(red, green, blue);
    elseif method == "linear"
        Ccam = bilinear(red, green, blue);
    else
        fprintf("Method not supported\n")
        Ccam =[];
    end
end


function Ccam = grbg(rawim, wbcoeffs, method, M, N)
    % % extract color components
    red = zeros(M, N);
    blue = zeros(M, N);
    red(1:2:end, 2:2:end) = rawim(1:2:end, 2:2:end);
    blue(2:2:end, 1:2:end) = rawim(2:2:end, 1:2:end);
    green = rawim - red - blue;

    wbmask = zeros(M, N);
    wbmask(1:2:end, 2:2:end) = wbcoeffs(1);
    wbmask(2:2:end, 1:2:end) = wbcoeffs(3);
    red = red .* wbmask;
    blue = blue .* wbmask;

    if method == "nearest"
        % grbg is essentially the same as gbrg with r and b flipped
        Ccam = nearestGBRG(blue, green, red);

        red = Ccam(:, :, 3);
        Ccam(:, :, 3) = Ccam(:, :, 1);
        Ccam(:, :, 1) = red;
    elseif method == "linear"
        Ccam = bilinear(red, green, blue);
    else
        fprintf("Method not supported\n")
        Ccam =[];
    end
end


function Ccam = nearestRGGB(red, green, blue)
    % MxN new size of image
    [M, N] = size(red);

    % % % fill gaps

    % % fill red on blue censors spots
    % nearest neighbor is considered NW value
    red(2:2:end, 2:2:end) = red(1:2:end, 1:2:end);

    % % fill red on green censors spots
    % nearest neighbor is considered W and N value
    red(1:2:end, 2:2:end) = red(1:2:end, 1:2:end);
    red(2:2:end, 1:2:end) = red(1:2:end, 1:2:end);

    % % fill green on red censors spots
    % nearest neighbor is considered E value
    green(1:2:end, 1:2:end) = green(1:2:end, 2:2:end);

    % % fill green on blue censors spots
    % nearest neighbor is considered W value
    green(2:2:end, 2:2:end) = green(2:2:end, 1:2:end);

    % % fill blue on red censors spots
    % nearest neighbor is considered SE value
    blue(1:2:end, 1:2:end) = blue(2:2:end, 2:2:end);

    % % fill blue on green censors spots
    % nearest neighbor is considered S and E value
    blue(1:2:end, 2:2:end) = blue(2:2:end, 2:2:end);
    blue(2:2:end, 1:2:end) = blue(2:2:end, 2:2:end);

    Ccam = zeros(M, N, 3);
    Ccam(:, :, 1) = red;
    Ccam(:, :, 2) = green;
    Ccam(:, :, 3) = blue;
end


function Ccam = nearestGBRG(red, green, blue)
    % MxN new size of image
    [M, N] = size(red);

    % % % fill gaps

    % % fill red on green censors spots
    % nearest neighbor is considered S and W value
    red(1:2:end, 1:2:end) = red(2:2:end, 1:2:end);
    red(2:2:end, 2:2:end) = red(2:2:end, 1:2:end);

    % % fill red on blue censors spots
    % nearest neighbor is considered SW value
    red(1:2:end, 2:2:end) = red(2:2:end, 1:2:end);

    % % fill green on red censors spots
    % nearest neighbor is considered E value
    green(2:2:end, 1:2:end) = green(2:2:end, 2:2:end);

    % % fill green on blue censors spots
    % nearest neighbor is considered W value
    green(1:2:end, 2:2:end) = green(1:2:end, 1:2:end);

    % % fill blue on red censors spots
    % nearest neighbor is considered NE value
    blue(2:2:end, 1:2:end) = blue(1:2:end, 2:2:end);

    % % fill blue on green censors spots
    % nearest neighbor is considered E and N value
    blue(1:2:end, 1:2:end) = blue(1:2:end, 2:2:end);
    blue(2:2:end, 2:2:end) = blue(1:2:end, 2:2:end);

    Ccam = zeros(M, N, 3);
    Ccam(:, :, 1) = red;
    Ccam(:, :, 2) = green;
    Ccam(:, :, 3) = blue;
end


function Ccam = bilinear(red, green, blue)
    % MxN new size of image
    [M, N] = size(red);

    % % fill gaps
    
    % fill red on blue spots
    redBlue = imfilter(red, [1 0 1; 0 0 0; 1 0 1;] / 4);

    % fill red on green spots
    redGreen1 = imfilter(red, [1 0 1;] / 2);
    redGreen2 = imfilter(red, [1; 0; 1] / 2);
    redGreen = redGreen1 + redGreen2;

    red = red + redBlue + redGreen;
    
    % green filling is the same for both red and blue
    greenRedBlue = imfilter(green, [0 1 0; 1 0 1; 0 1 0;] / 4);
    
    green = green + greenRedBlue;

    % fill blue on red spots
    blueRed = imfilter(blue, [1 0 1; 0 0 0; 1 0 1;] / 4);
    
    % fill blue on green spots
    blueGreen1 = imfilter(blue, [1 0 1;] / 2);
    blueGreen2 = imfilter(blue, [1; 0; 1] / 2);
    blueGreen = blueGreen1 + blueGreen2;

    blue = blue + blueRed + blueGreen;

    Ccam = zeros(M, N, 3);
    Ccam(:, :, 1) = red;
    Ccam(:, :, 2) = green;
    Ccam(:, :, 3) = blue;
end

