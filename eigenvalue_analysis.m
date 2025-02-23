% Parameters (Updated)
g0 = 4.2;          % Value for g_0
mu0 = 3.9;         % Value for mu_0
m = 3/4;           % Value for m
delta_t = 0.01;    % Time step (delta t)
t_max = 1;         % Maximum value for t
t_min = 0;         % Minimum value for t

% Compute w_values and w_count
w_min = 0.001;     % Minimum value for w
w_max = 100;       % Maximum value for w
w_count = log10(w_max / w_min) * 20 + 1;  % Number of points in w

% Generate w values logarithmically spaced
w_values = logspace(log10(w_min), log10(w_max), w_count); 

% Initialize the matrix M
M = zeros(w_count);

% Construct the matrix M with A_w^t and B_w^t terms
for w = 1:w_count
    % Compute A_w^t and B_w^t
    A_wt = 1 + g0 * w_values(w)^m * (delta_t / (w_max - w_min)) + mu0 * w_values(w)^(m-1) * delta_t;
    
    if w < w_count
        B_wt = -g0 * w_values(w)^m * (delta_t / (w_max - w_min)); % Only for w < w_max
    else
        B_wt = 0; % No B term for the last row
    end
    
    % Assign values to the matrix M
    M(w, w) = A_wt; % Diagonal elements are A_w^t
    if w < w_count
        M(w, w+1) = B_wt; % Superdiagonal elements are B_w^t
    end
end

% Perform eigenvalue analysis
[eigenvectors, eigenvalues] = eig(M);

% Extract eigenvalues and sort them
lambda = diag(eigenvalues);

% Display the eigenvalues
disp('Eigenvalues of the matrix M:');
disp(lambda);

% Plot the eigenvalues
figure;
plot(real(lambda), imag(lambda), 'bo', 'MarkerSize', 8);
xlabel('Real Part of Eigenvalues');
ylabel('Imaginary Part of Eigenvalues');
title('Eigenvalue Spectrum');
grid on;

% Check stability (magnitude of eigenvalues)
stable = all(abs(lambda) < 1); % Stable if all eigenvalues have magnitude < 1

if stable
    disp('The system is stable (all eigenvalues have magnitude < 1).');
else
    disp('The system is unstable (some eigenvalues have magnitude >= 1).');
end
