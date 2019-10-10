function [Ans, f, Aeq, beq, Res] = sudokuSolve(sudokuMat)
% This function can solve sudoku problems with any size, e.g., 4x4, 9x9.
% The size must be N^2 x N^2 (will be checked)
% Input: sodokuMat is the problem matrix. e.g.,
% [0 0 4 0;
%  2 0 0 1;
%  3 0 0 4;
%  0 2 0 0]
% Output: Ans is the completed sudoku matrix.
% f, Aeq, beq are params for intlinprog(). Res is the raw solution.

% ### The idea of the algorithm ###
% We construct a N^2 x N^2 x N^2 tensor x, where each element x(i, j, k) is
% binary (0 or 1). Each element x(i, j, m) in page m equals 1 if position
% (i, j) is digit m in sudoku matrix. Otherwise x(i, j, m) = 0.
% We can easily build constraints by following the sudoku rules.

M = sudokuMat';
% check if M is a square matrix
if size(M, 1) ~= size(M, 2)
    Ans = [];
    return
end

% check if M is N^2 by N^2
if sqrt(size(M, 1)) ~= int32(sqrt(size(M, 1)))
    Ans = [];
    return
end

N = sqrt(size(M, 1));


% constraints for known digits
valFixed = M(M > 0);
indFixed = find(M > 0);
numFixed = size(valFixed, 1);
I = floor(indFixed / N^2) + 1;
J = indFixed - (I - 1) * N^2;
A0 = zeros(numFixed, N^2 * N^2 * N^2);
ind = N^4 * (valFixed - 1) + (I - 1) * N^2 + J;
for t = 1 : numFixed
    A0(t, ind(t)) = 1;
end
b0 = ones(size(valFixed));

% page-wise constraints
A1 = zeros(N^2 * N^2, N^2 * N^2 * N^2);
t = 1;
for I = 1 : N^2
    for J = 1 : N^2
        ind = N^4 * ((1 : N^2) - 1) + (I - 1) * N^2 + J;
        A1(t, ind) = 1;
        t = t + 1;
    end
end

% constraints for columns
A2 = zeros(N^2 * N^2, N^2 * N^2 * N^2);
t = 1;
for K = 1 : N^2
    for J = 1 : N^2
        ind = N^4 * (K - 1) + ((1 : N^2) - 1) * N^2 + J;
        A2(t, ind) = 1;
        t = t + 1;
    end
end

% constraints for rows
A3 = zeros(N^2 * N^2, N^2 * N^2 * N^2);
t = 1;
for I = 1 : N^2
    for K = 1 : N^2
        ind = N^4 * (K - 1) + (I - 1) * N^2 + (1 : N^2);
        A3(t, ind) = 1;
        t = t + 1;
    end
end

% constraints for subgrids
A4 = zeros(N^2 * N^2, N^2 * N^2 * N^2);
t = 1;
for K = 1 : N^2
    for u = 0 : N : (N - 1) * N
        for v = 0 : N : (N - 1) * N
            for I = 1 : N
                for J = 1 : N
                    ind = N^4 * (K - 1) + (I + u - 1) * N^2 + J + v;
                    A4(t, ind) = 1;
                end
            end
            t = t + 1;
        end
    end
end

% solve Integer LP
f = zeros(N^2 * N^2 * N^2, 1);
Aeq = [A0; A1; A2; A3; A4];
beq = [b0; ones(size([A1; A2; A3; A4], 1), 1)];
lb = zeros(N^2 * N^2 * N^2, 1);
ub = ones(N^2 * N^2 * N^2, 1);
intcon = 1 : (N^2 * N^2 * N^2);
Res = intlinprog(f, intcon, [], [], Aeq, beq, lb, ub);

% interpret solutions
Ans = zeros(N^2, N^2);
for K = 1 : N^2
    for I = 1 : N^2
        for J = 1 : N^2
            ind = N^4 * (K - 1) + N^2 * (I - 1) + J;
            if Res(ind) > 0
                Ans(I, J) = K;
            end
        end
    end
end
return
end