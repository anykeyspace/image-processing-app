function [ ntsRatio ] = noiseToSignalRatio( imOriginal, imNoise )
%noiseToSignalRatio Определяет отношение шум-сигнал исходя из сравнения
%оригинальных данных с зашумленными.
%   Detailed explanation goes here

    noise = imNoise - imOriginal;
    mcolumn = mean(noise);
    noise = noise - mean(mcolumn);
    figure, imshow(noise);
    
    mcolumn = mean(imOriginal);
    data = imOriginal - mean(mcolumn);

    ntsRatio = mean(var(noise)) / mean(var(data));

end

