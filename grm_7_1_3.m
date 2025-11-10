% [[7,1,3:2]] Quantum RM code which can correct single X-type data and
% measurement errors and upto 3-data and measuremnt errors of Z type

[H_1, points] = kasami_grm_paritycheck_kasami(1, 3);
H_2 = H_1;
[m_1,n] = size(H_1);
m_2 = m_1;

k_1 = n-m_1;
k_2 = n-m_2;

%Ident_m = [eye(m_1), zeros(m_1, m_2);zeros(m_2, m_1), eye(m_2)];
% H_ds = [H_1, zeros(m_1,n), eye(m_1), zeros(m_1,m_2),zeros(m_1,r_1),zeros(m_1,r_2);
%     h1_per, zeros(r_1,n), zeros(r_1,m_1), zeros(r_1,m_2), eye(r_1), zeros(r_1,r_2);
%     zeros(m_2, n), H_2, zeros(m_2,m_1), eye(m_2), zeros(m_2,r_1), zeros(m_2, r_2);
%     zeros(r_2,n), h2_per, zeros(r_2, m_1), zeros(r_2, m_2), zeros(r_2, r_1), eye(r_2)
% ];



[h1_per, is_codewordx] = stack_left_shifts_in_code(H_1, 3,1); 
[h2_per, is_codewordz] = stack_left_shifts_in_code(H_2, 3,1); 
r_1 = size(h1_per,1);
r_2 = size(h2_per,1);

H_dsx = [H_1, eye(m_1), zeros(m_1,r_1);
          h1_per, zeros(r_1,m_1), eye(r_1)];

H_dsz = [H_2, eye(m_2), zeros(m_2,r_2);
          h1_per, zeros(r_2,m_2), eye(r_2)];


% File to save syndromes
filename = 'syndromes_7_1_3.mat';
save(filename, '-v7.3');
t1 = 1;
t2 = 1;
M1 = n+m_1+r_1;
M2 = n+m_2+r_2;
% Step 1: Generate and save syndromes incrementally
saveSyndromesToFile(M1, M2, t1, t2, H_dsx, H_dsz, filename);

% Step 2: Check for uniqueness of the saved syndromes
[countx, countz] = checkUniqueness(filename);

% Display results
if countx == 0
    disp('All X type syndromes are unique.');
else
    disp(['Found ' num2str(countx) ' duplicate X syndromes.']);
end

if countz == 0
    disp('All Z type syndromes are unique.');
else
    disp(['Found ' num2str(countz) ' duplicate Z syndromes.']);
end


% --------------------- Helper Functions ---------------------------
% Function to save syndromes to a .mat file incrementally
function saveSyndromesToFile(M1, M2, t1, t2, H_dsx, H_dsz, filename)
    % Ensure the file exists before appending
    if ~isfile(filename)
        save(filename, '-v7.3');
    end

    % Generate and save bit-flip syndromes in parallel
    parfor i = 1:t2
        % Generate bit-flip error patterns of weight i
        bitFlipPatterns = generateErrorPatterns(M2, i);
        
        % Calculate syndromes for these bit-flip error patterns
        bitFlipSyndromes = mod(H_dsz * bitFlipPatterns', 2);

        % Save syndromes locally for each worker
        saveBitFlipSyndromes(bitFlipSyndromes, i, filename);
    end

    % Generate and save phase-flip syndromes in parallel
    parfor i = 1:t1
        % Generate phase-flip error patterns of weight i
        phaseFlipPatterns = generateErrorPatterns(M1, i);
        
        % Calculate syndromes for these phase-flip error patterns
        phaseFlipSyndromes = mod(H_dsx * phaseFlipPatterns', 2);

        % Save syndromes locally for each worker
        savePhaseFlipSyndromes(phaseFlipSyndromes, i, filename);
    end

    disp('Syndrome generation and saving completed.');
    combineSyndromes(filename, t1,t2);
end

function combineSyndromes(filename, t1, t2)
    % Access the MAT-file for appending
    file = matfile(filename, 'Writable', true);

    % Initialize or retrieve combined bit-flip syndromes
    if ~isprop(file, 'bitFlipSyndromes')
        file.bitFlipSyndromes = [];  % Initialize if not already present
    end

    % Combine all bit-flip syndromes
    for i = 1:t2
        % Load individual syndromes dynamically
        varName = ['bitFlipSyndromes_weight_', num2str(i)];
        tempData = load(filename, varName);
        
        % Append the loaded syndromes
        if isfield(tempData, varName)
            file.bitFlipSyndromes = [file.bitFlipSyndromes, tempData.(varName)];
        else
            disp(['Variable ', varName, ' not found in file.']);
        end
    end
    disp('Combined all bit-flip syndromes.');

    % Initialize or retrieve combined phase-flip syndromes
    if ~isprop(file, 'phaseFlipSyndromes')
        file.phaseFlipSyndromes = [];  % Initialize if not already present
    end

    % Combine all phase-flip syndromes
    for i = 1:t1
        % Load individual syndromes dynamically
        varName = ['phaseFlipSyndromes_weight_', num2str(i)];
        tempData = load(filename, varName);
        
        % Append the loaded syndromes
        if isfield(tempData, varName)
            file.phaseFlipSyndromes = [file.phaseFlipSyndromes, tempData.(varName)];
        else
            disp(['Variable ', varName, ' not found in file.']);
        end
    end
    disp('Combined all phase-flip syndromes.');
end

% Function to save bit-flip syndromes to a file (worker-specific chunks)
function saveBitFlipSyndromes(bitFlipSyndromes, weight, filename)
    varName = ['bitFlipSyndromes_weight_', num2str(weight)];
    S.(varName) = bitFlipSyndromes; % Create a struct for dynamic variable naming
    save(filename, '-struct', 'S', '-append');
    disp(['Saved bit-flip syndromes for weight ', num2str(weight)]);
end

% Helper function for saving phase-flip syndromes
function savePhaseFlipSyndromes(phaseFlipSyndromes, weight, filename)
    varName = ['phaseFlipSyndromes_weight_', num2str(weight)];
    S.(varName) = phaseFlipSyndromes; % Create a struct for dynamic variable naming
    save(filename, '-struct', 'S', '-append');
    disp(['Saved phase-flip syndromes for weight ', num2str(weight)]);
end


% Function to generate all error patterns of length n with specified Hamming weight
function patterns = generateErrorPatterns(n, weight)
    indices = sparse(nchoosek(1:n, weight));  % All combinations of positions with '1's
    patterns = sparse(zeros(size(indices, 1), n));
    for i = 1:size(indices, 1)
        patterns(i, indices(i, :)) = 1;  % Set the positions as '1'
    end
end

% Function to check uniqueness of syndromes from the saved .mat file
function [countx, countz] = checkUniqueness(filename)
    % Initialize counts for non-unique syndromes
    countx = 0;
    countz = 0;

    % Access the file using matfile for efficient partial loading
    data = matfile(filename);

    % Check uniqueness of bit-flip syndromes
    if isprop(data, 'bitFlipSyndromes')
        bitFlipSyndromes = data.bitFlipSyndromes;
        % Use unique to check for duplicate columns
        [~, ia, ~] = unique(bitFlipSyndromes', 'rows', 'stable');
        countx = size(bitFlipSyndromes, 2) - numel(ia);  % Non-unique count
    else
        disp('bitFlipSyndromes not found in file.');
    end

    % Check uniqueness of phase-flip syndromes
    if isprop(data, 'phaseFlipSyndromes')
        phaseFlipSyndromes = data.phaseFlipSyndromes;
        % Use unique to check for duplicate columns
        [~, ia, ~] = unique(phaseFlipSyndromes', 'rows', 'stable');
        countz = size(phaseFlipSyndromes, 2) - numel(ia);  % Non-unique count
    else
        disp('phaseFlipSyndromes not found in file.');
    end
end