function C = codewords_from_G(G)
    % G: Generator matrix (k x n)
    k = size(G, 1);          % Number of rows of G = length of input bitwords
    all_msgs = de2bi(0:2^k-1, k, 'left-msb'); % All 2^k input vectors
    C = mod(all_msgs * G, 2); % Multiply and reduce mod 2 to get codewords
end
