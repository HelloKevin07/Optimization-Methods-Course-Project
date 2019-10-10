% We use a 4x1 binary vector Pt to represent the positions of alpha, beta, 
% gamma, delta at time t. We use 0 to represent the left side of the river, 
% 1 for the right side, e.g., we have P0 = [0 0 0 0] and P5 = [1 1 1 1].
%
% So the LP problem is
% ( 1' means [1 1 1 1] transpose, T = [1 2 5 10])
%
%  min   c1 + c2 + c3 + c4 + c5
%  s.t.   1'P1 = 2
%        -1'(P2 - P1) = 1
%         1'(P3 - P2) = 2
%        -1'(P4 - P3) = 1
%         1'(P5 - P4) = 2
%
%         T.* P1 <= c1
%        -T.*(P2 - P1) <= c2
%         T.*(P3 - P2) <= c3
%        -T.*(P4 - P3) <= c4
%         T.*(P5 - P4) <= c5
%         
%        -P1 <= 0
%         P2 - P1 <= 0
%         P2 - P3 <= 0
%         P4 - P3 <= 0
%         P4 - P5 <= 0
%         
%         c1, c2, c3, c4, c5, P1, P2, P3, P4, P5 >= 0
%         P1, P2, P3, P4, P5 <= 1
%
% There are 25 unkonwn variables, i.e., x = [P1; P2; P3; P4; P5; c1; c2; c3; c4; c5].
% To solve this LP, we must use intlinprog(), otherwise we won't get 
% integer solutions. The optimal solution is explained at the end. Please 
% roll down and run the code.

f = [zeros(20, 1); ones(5, 1)];
% first part of constraints
Aeq1 = [ones(4, 1); zeros(21, 1);];
Aeq2 = -[-ones(4, 1); ones(4, 1); zeros(17, 1)];
Aeq3 = [zeros(4, 1); -ones(4, 1); ones(4, 1); zeros(13, 1)];
Aeq4 = -[zeros(8, 1); -ones(4, 1); ones(4, 1); zeros(9, 1)];
Aeq5 = [zeros(12, 1); -ones(4, 1); ones(4, 1); zeros(5, 1)];
Aeq = [Aeq1'; Aeq2'; Aeq3'; Aeq4'; Aeq5'];
beq = [2 1 2 1 2]';

% second part of constraints
c = [1 2 5 10];
A1 = zeros(20, 25);
for k = 1 : 5
    for u = 1 : 4
        A1((k - 1) * 4 + u, 20 + k) = -1;
        A1((k - 1) * 4 + u, (k - 1) * 4 + u) = c(u) * (-1)^(k - 1);
        if k ~= 1
            A1((k - 1) * 4 + u, (k - 1) * 4 + u - 4) = -c(u) * (-1)^(k - 1);
        end
    end
end

% third part of constraints
A2 = zeros(20, 25);
for k = 1 : 5
    for u = 1 : 4
        A2((k - 1) * 4 + u, (k - 1) * 4 + u) = -(-1)^(k - 1);
        if k ~= 1
           A2((k - 1) * 4 + u, (k - 1) * 4 + u - 4) = (-1)^(k - 1); 
        end
    end
end

A = [A1; A2];

b = zeros(40, 1);

% fourth part of constraints
lb = zeros(25, 1);
ub = [ones(20, 1); Inf * ones(5, 1)];
intcon = 1 : 25;

% solve the integer LP
x = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);

% interpret the solution
P1 = x(1 : 4)';
P2 = x(5 : 8)';
P3 = x(9 : 12)';
P4 = x(13 : 16)';
P5 = x(17 : 20)';

% position matrix P
P = int32([P1; P2; P3; P4; P5]);

S1 = P1;
S2 = P2 - P1;
S3 = P3 - P2;
S4 = P4 - P3;
S5 = P5 - P4;

% transition matrix S
S = int32([S1; S2; S3; S4; S5]);

fprintf('The position matrix P is\n')
disp(P)
fprintf('If we represent each transition at time t by St = Pt - Pt-1\n')
fprintf('Then the transition matrix S is\n')
disp(S)
% From the transition matrix S, we can know the optimal scheme is 
% (1) alpha + beta cross foward  --- 2 min
% (2) beta cross back            --- 2 min
% (3) gamma + delta cross foward --- 10 min
% (4) alpha cross back           --- 1 min
% (5) alpha + beta cross foward  --- 2 min
% Total time = 17 min
 





