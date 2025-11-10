function out = is_codeword_of_G(v, G)
    % v: 1 x n binary row vector
    % G: k x n binary generator matrix
    v = mod(v,2); G = mod(G,2);

    % Try to solve G' * x = v' over GF(2)
    x = gflineq(G', v', 2); % Returns a k x 1 solution or NaN if no solution exists
    out = ~any(isnan(x));   % True if solution exists, false if not
end
