function sudokuGrid = generate4x4Sudoku()
    % Define the basic solution for a 4x4 Sudoku grid
    baseGrid = [1, 2, 3, 4;
                3, 4, 1, 2;
                2, 1, 4, 3;
                4, 3, 2, 1];
            
    % Randomly permute rows and columns within blocks
    perm = randperm(2); % For the 2x2 blocks
    rowPerm = [perm, perm + 2];
    colPerm = [perm, perm + 2];
    
    % Apply row and column permutations
    sudokuGrid = baseGrid(rowPerm, colPerm);
    
    % Remove a few numbers randomly to create a puzzle
    numRemove = 6; % Number of cells to clear
    for i = 1:numRemove
        row = randi([1, 4]);
        col = randi([1, 4]);
        sudokuGrid(row, col) = 0; % Use 0 to indicate empty cells
    end
    
    % Display the puzzle grid
    %disp('4x4 Sudoku Puzzle (0 indicates an empty cell):');
    %disp(sudokuGrid);
end

% Generate and display the Sudoku puzzle
%sudokuGrid = generate4x4Sudoku();
