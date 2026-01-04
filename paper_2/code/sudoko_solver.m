clc;
clear;
%sudoku problem 
sudoku_matrix =generate9x9Sudoku();
%vectorize the matrix
sudoku_matrix_tp=transpose(sudoku_matrix);
sudoku_vector=sudoku_matrix_tp(:);
%% constraint setup
% Define the size of each block
N=size(sudoku_matrix,2);
% Create the identity and zero blocks
I_block = eye(N);
J_block=repmat(I_block,1,sqrt(N));
zero_block = zeros(N);
%% row compatibility
rows = N^2;
cols = N^3;

% Initialize the matrix with zeros
matrix_row = zeros(rows, cols);

% Fill in the blocks as specified
for i = 0:N-1  % for each 9-row block
    start_row = i * N + 1;
    for j = 0:N-1  % place 9 consecutive identity matrices in each block
        start_col = i * N^2 + j * N + 1;
        matrix_row(start_row:start_row + N-1, start_col:start_col + N-1) = eye(N);
    end
end
%% Column compatibility
% Initialize a cell array to store the matrices
matrices_col = cell(1, N);

% Generate each matrix with sliding eye(N)
for i = 1:N
    % Initialize a 9x81 matrix with zeros
    M = zeros(N, N^2);
    
    % Define the starting column for the eye(9) block in the current matrix
    start_col = (i - 1) * N + 1;
    
    % Place eye(9) in the designated column range
    M(:, start_col:start_col + N-1) = eye(N);
    
    % Store the matrix in the cell array
    matrices_col{i} = M;
end
matrix_column=[];
for i=1:N
matrix_col=repmat(matrices_col{i},1,N);
matrix_column=[matrix_column;matrix_col];
end
%% constraint from box
% Calculate sqrt(N)
sqrt_N = sqrt(N);

% Define identity and zero matrices
I_N = eye(N);                    % NxN identity matrix
J_N = repmat(I_N, 1, sqrt_N);    % Nx(sqrt(N)*N) matrix with sqrt(N) identity matrices side by side
Z_N= zeros(N);
Z_2N = repmat(Z_N,1,sqrt_N); % Nx(N*(sqrt(N)-1)) zero matrix
% Construct the repeating block: [J_N, Z_2N]
block = [J_N, Z_2N];  % This block has size N x N^2
% Initialize a cell array to store sqrt(N) x sqrt(N) blocks
matrix_box = [];
% Fill the cell matrix with the appropriate blocks
for i = 1:sqrt_N   
    for j = 1:sqrt_N
         current_block = zeros(N, N^3);
        % Create a new block with J_N positioned correctly
        start_col1 = (i - 1) * N*size(Z_2N,2)+(j-1)*size(Z_2N,2)+1; % Position for J_N within the block
        for k=1:sqrt_N
        start_col=start_col1+(k-1)*(size(Z_2N,2)*(sqrt_N-1)+size(J_N,2));
        current_block(:, start_col:start_col + size(J_N, 2) - 1) = J_N;
        end
        matrix_box=[matrix_box;current_block]; 
    end
    
end
%% constraint from clue
%comb=[J_block,repmat(zero_block,1,N-sqrt(N))];
index=find(sudoku_vector); 
matrix_clue=zeros(N^2,N^3);
u1=zeros(1,N);
u2=eye(N);
u=[u1;u2];
for m=1:N^2
    x=sudoku_vector(m);
    startIdx=(m-1)*N+1;
    endIdx=m*N;
    matrix_clue(m,:) = [matrix_clue(m,1:startIdx-1), u(x+1,:), matrix_clue(m,endIdx+1:end)];
end
v=matrix_clue;
matrix_clue(all(matrix_clue == 0, 2), :) = [];
%% constraint from blank
matrix_cell=zeros(N^2,N^3);
for m=1:N^2
    matrix_cell(m,(m-1)*N+1:(m-1)*N+N)=1;
end
%% combined_constraint
A=[matrix_row;matrix_column;matrix_box;matrix_cell;matrix_clue];
b = ones(size(A, 1), 1); % Vector of ones with the same number of rows as A
%% optimization_problem
cvx_begin
variable y(size(A, 2), 1);
    minimize(norm(y, 1))   % Minimize the L1 norm of x
    subject to
        A * y == b         % Constraint Ax = b
cvx_end
%% get the solved sudoku
x=round(y);
sudoku_solved=zeros(N^2,N);
for i=1:N^2
    for j=1:N
    sudoku_solved(i,j)=x(N*(i-1)+j);
    end
end
%%
sudoku_vector_solved=zeros(N^2,1);
for i=1:N^2
    sudoku_vector_solved(i)=find(sudoku_solved(i,:));
end
%%
sudoku_matrix_solved=zeros(N,N);
for i=1:N
    for j=1:N
        sudoku_matrix_solved(i,j)= sudoku_vector_solved((N*(i-1)+j));
    end
end