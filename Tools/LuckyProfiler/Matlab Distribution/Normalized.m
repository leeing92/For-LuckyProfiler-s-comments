function [oI] = Normalized(I)
%this function is used to normalized the image for scale range and type
%   此处显示详细说明
    minI=min(min(I));
    maxI=max(max(I));
    oI=(I-minI)/(maxI-minI);
% oI=I/max(max(I));
end

