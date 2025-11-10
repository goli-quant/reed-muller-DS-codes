function A = findA(H_prime, H)
    % Find binary matrix A with H_prime = A * H (mod 2)
    % H_prime: t x n, H: k x n
    % Output: A is t x k binary

    [t, n] = size(H_prime);
    k = size(H, 1);
    A = zeros(t, k);

    for i = 1:t
        % For each row of H_prime, find binary vector alpha with alpha*H = H_prime(i,:) mod 2
        alpha = solve_gf2_linear(H, H_prime(i,:));
        A(i, :) = alpha;
    end
end

function alpha = solve_gf2_linear(H, v)
    % For small k, brute-force all possible binary vectors alpha (length k)
    k = size(H,1);
    all_alpha = de2bi(0:2^k-1, k, 'left-msb');
    for idx = 1:size(all_alpha,1)
        if isequal(mod(all_alpha(idx,:) * H, 2), v)
            alpha = all_alpha(idx,:);
            return;
        end
    end
    error('No binary solution for row!');
end
