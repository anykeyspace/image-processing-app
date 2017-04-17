function [ k ] = noiseToSignalRatio2( imOriginal, imNoise, imFiltered )
%noiseToSignalRatio2
%   Detailed explanation goes here

    a = std2(imNoise - imOriginal);

    b = std2(imFiltered - imOriginal);
    
    k = a / b;

end

