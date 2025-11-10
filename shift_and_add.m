function [H_prime, is_codeword] = shift_and_add(H, m, r)
% H: k x n (parity-check matrix, k = size(H,1), n = size(H,2))
% m: # of degree-1 generators (first m rows are cyclic block)
% r: code order (not used directly, included for clarity)

[p,n] = size(H);

all_ones_vec = ones(1, n);

H_prime = zeros(p, n);
is_codeword = false(1, p);

first = H(1,:);
shifted = first;
for i = 1:m-1
    shifted = circshift(shifted, [0, -1]);
end
% Block 1: Degree-1 (first m rows), m cyclic shifts (each row is a circshift)
for i=1:m
    shifted = circshift(shifted, [0, -1]);
    H_prime(i, :) = mod(shifted + all_ones_vec, 2);
    is_codeword(i) = is_codeword_of_G(shifted, H);
end

temp = m;
if r>1
    for i=2:r
        current_block_size = nchoosek(m,i);
        current_pointer = temp+current_block_size;
        shifted = H(current_pointer-1,:);
        % Remaining blocks: higher-degree (rows m+1 to k), just cyclic shift own rows
        for j = 1:current_block_size
            shifted = circshift(shifted, [0, -1]); 
            H_prime(temp, :) = shifted;
            is_codeword(j) = is_codeword_of_G(shifted, H);
        end
    end
end
end
