function [ outImage ] = SVDImageProcessing_impulse( inImage )
% SVDImageProcessing_impulse
%   [ outImage ] = SVDImageProcessing_impulse( inImage )
    
%    k = 3; % Количество элементов в передаваемой в фильтр строке

    [m, n, p] = size(inImage);
    outImage = zeros(m, n, p);
    
%     for rgb = 1 : p
%         for i = 1 : m
%             for j = 1 : k : n
%                 if n - j - k > 0
%                     outImage(i, j : j+k, rgb) = ImpulseSVDfilter( inImage(i, j : j+k, rgb) );
%                 else
%                     outImage(i, j : n, rgb) = ImpulseSVDfilter( inImage(i, j : n, rgb) );
%                 end
%             end
%         end
%     end
    
    for rgb = 1 : p
        for i = 1 : n
            outImage(:, i, rgb) = ImpulseSVDfilter( inImage(:, i, rgb) )';
        end
    end

end

function [ oStr ] = ImpulseSVDfilter( iStr )
% ImpulseSVDfilter
%   [ oStr ] = ImpulseSVDfilter( iStr )

	n = length(iStr);
    X(1, :) = iStr(1 : n-1);
    X(2, :) = iStr(2 : n);
    
    [U, S, V] = svd( X );
%    a1 = U(:, 1) * S(1, 1) * V(:, 1)';
    a2 = U(:, 2) * S(2, 2) * V(:, 2)';
	
	noise = cat(2, a2(1, :), a2(2, n - 1));
    
	% high and low values of noise
	nh = max(iStr);
	nl = min(iStr);
	
	% removing noise from the input signal
	oStr = zeros(1, n);
	for i = 1 : n
		if iStr(i) == nh || iStr(i) == nl
			oStr(i) = iStr(i) - 2 * noise(i);
		else
			oStr(i) = iStr(i);
		end
	end
    
end

