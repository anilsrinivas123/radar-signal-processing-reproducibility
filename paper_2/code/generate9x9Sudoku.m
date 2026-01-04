function sudokuGrid = generate9x9Sudoku()
    % Define the basic solution for a 9x9 Sudoku grid
    baseGrid = [
        1 2 3  4 5 6  7 8 9;
        4 5 6  7 8 9  1 2 3;
        7 8 9  1 2 3  4 5 6;
        
        2 3 4  5 6 7  8 9 1;
        5 6 7  8 9 1  2 3 4;
        8 9 1  2 3 4  5 6 7;
        
        3 4 5  6 7 8  9 1 2;
        6 7 8  9 1 2  3 4 5;
        9 1 2  3 4 5  6 7 8
    ];

    % Permute rows and columns within each 3x3 block
    rowPerm = [randperm(3), randperm(3) + 3, randperm(3) + 6];
    colPerm = [randperm(3), randperm(3) + 3, randperm(3) + 6];
    
    % Apply row and column permutations
    sudokuGrid = baseGrid(rowPerm, colPerm);
    
    % Remove a few numbers randomly to create a puzzle
    numRemove = 40; % Number of cells to clear for puzzle difficulty
    for i = 1:numRemove
        row = randi([1, 9]);
        col = randi([1, 9]);
        sudokuGrid(row, col) = 0; % Use 0 to indicate empty cells
    end
end
