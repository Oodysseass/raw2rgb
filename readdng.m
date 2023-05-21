function [rawim, XYZ2Cam, wbcoeffs] = readdng(filename)
    % open .dng
    obj = Tiff(filename, 'r');
    offsets = getTag(obj, 'SubIFD');
    setSubDirectory(obj, offsets(1))
    rawim = read(obj);
    close(obj);
    
    % % get info
    meta_info = imfinfo(filename);
    
    % size info
    y_origin = meta_info.SubIFDs{1}.ActiveArea(1) + 1;
    x_origin = meta_info.SubIFDs{1}.ActiveArea(2) + 1;
    
    width = meta_info.SubIFDs{1}.DefaultCropSize(1);
    height = meta_info.SubIFDs{1}.DefaultCropSize(2);

    % black-white info
    blacklevel = meta_info.SubIFDs{1}.BlackLevel(1);
    whitelevel = meta_info.SubIFDs{1}.WhiteLevel;

    % white balance info
    wbcoeffs = (meta_info.AsShotNeutral) .^ -1;
    wbcoeffs = wbcoeffs / wbcoeffs(2);
    
    % color space info
    XYZ2Cam = meta_info.ColorMatrix2;
    XYZ2Cam = reshape(XYZ2Cam, 3, 3)';

    % % correct info
    % keep true size
    rawim = double(rawim(y_origin : y_origin + height - 1, ...
                    x_origin : x_origin + width - 1));

    % subtract possible black offset due to noise
    rawim = rawim - blacklevel;
    % normalize 0-1
    rawim = rawim / (whitelevel - blacklevel);
    % correct out of bounds  values
    rawim = max(0, min(1, rawim));

    % % keep size as an even number
    if mod(2, size(rawim, 1)) == 1
        rawim(end, :) = [];
    end
    if mod(2, size(rawim, 2)) == 1
        rawim(:, end) = [];
    end
end

