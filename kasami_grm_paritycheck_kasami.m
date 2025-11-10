function [H, points] = kasami_grm_paritycheck_kasami(r, m)
    % Construct parity-check matrix H as per Kasami for cyclic RM(r,m)
    % Excludes all-ones vector; uses only monomials of degree >= 1 and <= m-1-r

    n = 2^m - 1; % Cyclic code length
    prim_poly = primpoly(m, 'min');
    alpha = gf(2, m, prim_poly);
    elements = gf(zeros(n,1), m, prim_poly);
    elements(1) = gf(1, m, prim_poly);
    for i = 2:n
        elements(i) = alpha^(i-1);  % Powers of alpha (nonzero field points)
    end
    vals = elements.x;
    A = de2bi(vals, m, 'left-msb')'; % m x n array: binary basis expansions
    points = A';

    r_dual = m - 1 - r;
    mon_idx_dual = {};
    for deg = 1:r_dual             % Only degrees 1 through r_dual
        combos = nchoosek(1:m, deg);
        for idx = 1:size(combos, 1)
            mon_idx_dual{end+1} = combos(idx, :);
        end
    end

    H = zeros(length(mon_idx_dual), n);
    for i = 1:length(mon_idx_dual)
        H(i,:) = prod(A(mon_idx_dual{i},:), 1);  % Evaluate selected monomials only
    end
    H = mod(H, 2);              % Modulo 2 for binary code
    % The resulting H is always exactly n-k in size; all rows are independent and exclude the all-ones vector.
end
