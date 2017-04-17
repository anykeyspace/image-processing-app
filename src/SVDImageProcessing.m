function [ outImage ] = SVDImageProcessing( inImage, maskSize, level )
% SVDImageProcessing
%   [ outImage ] = SVDImageProcessing( inImage, maskSize, level )

    [m, n, p] = size(inImage);
    outImage = inImage;
    halfMaskSize = floor(maskSize / 2);
    
    for rgb = 1 : p
        for x = 1 + halfMaskSize : m - halfMaskSize
            for y = 1 + halfMaskSize : n - halfMaskSize
                tempImage = inImage(x - halfMaskSize : x + halfMaskSize, y - halfMaskSize : y + halfMaskSize, rgb);
                tempImage = SVDfilter(tempImage, level);
                outImage(x, y, rgb) = tempImage(halfMaskSize + 1, halfMaskSize + 1);
            end
        end
    end

end

function [ outImage ] = SVDfilter( inImage, level )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    [m, n] = size(inImage);
    N = min(m, n);
    [U, S, V] = svd(inImage);

    rank = rankSVD(inImage, level);

    tempImage = zeros(m, n);
    
    S1 = log(1 + double(S));

    for i = 1 : rank
        tempImage = tempImage + U(:,i) * S1(i,i) * V(:,i)';
    end

    outImage = tempImage;

end

function [ rank ] = rankSVD( A, level )
% myRank возвращает ранг аппроксимирующей матрицы
%

    s = svd(A);
    
    Lopt = level; % Порог аппроксимации
    A_approx = s(1)^2;
    A_norm = sum( s.^2 );
    Lp = A_approx / A_norm;
    
    i = 2;
    while Lp < Lopt
        A_approx = A_approx + s(i)^2;
        Lp = A_approx / A_norm;
        i = i + 1;
    end
    
    rank = i - 1;
    
end
