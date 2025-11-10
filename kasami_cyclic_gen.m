m = 4;
n = 2^m - 1;
prim_poly_decimal = primpoly(m, 'min');
alpha = gf(2, m, prim_poly_decimal);

elements = gf(zeros(n,1), m, prim_poly_decimal);
for i = 1:n
    elements(i) = alpha^i;
end

G = ones(m+1, n);  % Generator matrix; first row is constant monomial
for i = 1:n
    val = elements(i).x;                           % Get scalar integer
    binvec = de2bi(val, m, 'left-msb');            % Get basis coefficients as bits
    G(2:end, i) = binvec.';                        % Stack as columns in G
end

disp('Generator matrix G:')
disp(G)
