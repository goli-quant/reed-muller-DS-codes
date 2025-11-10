% r = 1;
% m = 4;
% 
% 
% [G, points] = kasami_grm_matrix(r, m);
% [H, points] = kasami_grm_paritycheck_kasami(r, m);
% 
% [H_prime, is_codeWard] = shift_and_add_ones(H, m, m-r-1);

H = [0 0 1 0 1 1 1; 0 1 0 1 1 1 0; 1 0 0 1 0 1 1];
H_prime = [0 0 1 0 1 1 1; 0 1 0 1 1 1 0; 1 0 0 1 0 1 1; 1 1 1 0 0 1 0; 0 1 1 1 0 0 1];

A = findA(H_prime, H);