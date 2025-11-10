function isCyclicLeft = check_left_cyclicity(G)
    [k, n] = size(G);
    numWords = 2^k;
    codewords = zeros(numWords, n);
    for idx = 0:numWords-1
        v = de2bi(idx, k, 'left-msb')';
        codewords(idx+1, :) = mod(v' * G, 2);
    end

    isCyclicLeft = true;
    for i = 1:numWords
        cw = codewords(i, :);
        cw_shift = [cw(2:end), cw(1)];   % Left cyclic shift
        % Check if shifted codeword is in codewords
        found = false;
        for j = 1:numWords
            if isequal(cw_shift, codewords(j, :))
                found = true;
                break;
            end
        end
        if ~found
            isCyclicLeft = false;
            return;
        end
    end
end
