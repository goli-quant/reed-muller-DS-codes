function min_distance = greedy_min_distance(G_s)
    % Input:
    % G_s: Generator matrix of the shortened RM code (binary matrix)
    % Output:
    % min_distance: Estimated minimum distance of the code

    [k, n] = size(G_s);   % k = number of rows, n = length of the code
    min_distance = n;     % Initialize minimum distance as maximum possible

    % Iterate over all rows to start the greedy search
    for i = 1:k
        % Start with the current row as the initial codeword
        current_codeword = G_s(i, :);
        current_weight = sum(current_codeword); % Hamming weight of codeword

        % Update minimum distance if a new lower weight is found
        if current_weight > 0 && current_weight < min_distance
            min_distance = current_weight;
        end

        % Combine current row with all other rows greedily
        for j = i+1:k
            % Form a new codeword by combining two rows (binary addition)
            new_codeword = mod(current_codeword + G_s(j, :), 2);
            new_weight = sum(new_codeword); % Compute Hamming weight

            % Update minimum distance if a new lower weight is found
            if new_weight > 0 && new_weight < min_distance
                min_distance = new_weight;
            end

            % Continue combining further rows (nested greedy search)
            for l = j+1:k
                next_codeword = mod(new_codeword + G_s(l, :), 2);
                next_weight = sum(next_codeword);

                if next_weight > 0 && next_weight < min_distance
                    min_distance = next_weight;
                end
            end
        end
    end

    % Output the estimated minimum distance
    fprintf('Estimated minimum distance: %d\n', min_distance);
end