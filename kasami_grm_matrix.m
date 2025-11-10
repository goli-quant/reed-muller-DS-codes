function [G, points] = kasami_grm_matrix(r, m)
    % Outputs generator G and parity-check H matrices for cyclic RM(r,m) code (binary)
    % G: generator matrix (k x n)
    % H: parity-check matrix (n-k x n)
    % points: binary GF(2^m) vectors for each coordinate

    n = 2^m - 1;                                 % Cyclic code length
    prim_poly = primpoly(m, 'min');               % Primitive polynomial for GF(2^m)
    alpha = gf(2, m, prim_poly);
    elements = gf(zeros(n,1), m, prim_poly);      % List nonzero elements of GF(2^m)
    elements(1) = gf(1, m, prim_poly);
    for i = 2:n, elements(i) = alpha^(i-1); end
    vals = elements.x;
    A = de2bi(vals, m, 'left-msb')';              % m x n array, basis expansion
    points = A';

    % --- Generator matrix construction (G) ----
    mon_idx = {};                     % Indices of variable sets for monomials deg <= r
    for deg = 0:r
        cmb = nchoosek(1:m, deg);
        for i = 1:size(cmb,1)
            mon_idx{end+1} = cmb(i,:);
        end
    end
    k = length(mon_idx);
    G = zeros(k, n);
    for i = 1:k
        if isempty(mon_idx{i})
            G(i,:) = ones(1, n);
        else
            G(i,:) = prod(A(mon_idx{i},:), 1);
        end
    end
    G = mod(G, 2);
end
