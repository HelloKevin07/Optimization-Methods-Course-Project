% Suppose we have the following sudoku problem, where each '0' needs to be
% filled with 1~N digits.
Problem = ...
    [0 0 4 0;
     2 0 0 1;
     3 0 0 4;
     0 2 0 0];
% solve sudoku using Integer LP
[sudokuAns, c, Aeq, beq, res] = sudokuSolve(Problem);
% test the uniqueness of the solution
isuniq = uniquenessTest(c, Aeq, beq, res);
% show the compeleted sudoku solution
disp(sudokuAns);
% show the uniqueness
if isuniq
    fprintf('The answer is unique.\n');
else
    fprintf('The answer is NOT unique.\n');
end