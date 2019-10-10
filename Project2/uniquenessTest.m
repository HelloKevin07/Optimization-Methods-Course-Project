function [res] = uniquenessTest(c, A, b, xp)
% This function can test uniqueness of the soduku solution.
% Input: c, A, b is the params for integer LP of sudoku
% xp is the raw solution of the sudoku
% Output: res will be 1 if the solution is unique, else be 0.

% ### The idea of the algorithm ###
% To test the uniqueness of basic optimal solution xp of the following LP:
%     min cTx
%     s.t. Ax = b
%          x >= 0
% We need to solve the following LP(*):
%     max dTx
%     s.t. Ax = b
%          cTx = cTxp
%          x >= 0
% where d(j) = 1 for xp(j) = 0, d(j) = 0 for xp(j) > 0.
% If the optimum of the LP(*) is 0, we can conclude the solution is
% unique, otherwise it's not unique.
%
% The intuition for this algorithm is, we create a new LP, the feasible set
% of which contains the xp and other possible optimal solutions. However, 
% the objective we designed is to find a different basic feasible solution.
% If there is a differen basic optimal solution, it will cause dTx > 0.

d = (xp == 0);
Aeq = [A; c'];
beq = [b; c' * xp];
lb = zeros(size(xp));
ub = ones(size(xp));
Ans = linprog(-d, [], [], Aeq, beq, lb, ub);
res = (d' * Ans) == 0;
return
end