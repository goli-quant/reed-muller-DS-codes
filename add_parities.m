function H_rm = add_parities(H) 
  
    parities = mod(sum(H, 2), 2);  % Column vector of row parities

    H_rm = [parities, H];          % Concatenate to the left
end
