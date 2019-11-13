syms w a0 a1 a2 alpha0 alpha1 alpha2 
N = 5;
w = [0.5236; 1.0472; 1.5708; 2.0944; 2.6180];
Pw = [0.6821; 4.3232; 3.7540; 0.4368; 0.1988];
E(a0, a1, a2) = sym(0);

for i = 1 : 5
    Pwi_hat(a0, a1, a2) = (a0 + a1 * exp(-1j * w(i)) + a2 * exp(-2j * w(i)))...
        * conj(a0 + a1 * exp(-1j * w(i)) + a2 * exp(-2j * w(i)));
    tmp = Pw(i) * Pwi_hat - log(Pw(i) * Pwi_hat) - 1;
    E(a0, a1, a2) = E(a0, a1, a2) + tmp;
end

step = Inf;
distance = Inf;

a = [1 0 0]';

dE_a0 = vpa(diff(E, a0), 10);
dE_a1 = vpa(diff(E, a1), 10);
dE_a2 = vpa(diff(E, a2), 10);
count = 0;
v = 0;

tau = 0.9;
beta = 0.5;
lr = 1;

E_record = [];
distance_record = [];
optimal = [1; -0.5161; 0.9940];
count = 0;
while step > 1e-10
    
    E_pre = double(E(a(1), a(2), a(3)));
    
    grad_a0 = double(dE_a0(a(1), a(2), a(3)));
    grad_a1 = double(dE_a1(a(1), a(2), a(3)));
    grad_a2 = double(dE_a2(a(1), a(2), a(3)));
    
    % line search
    grad = [grad_a0; grad_a1; grad_a2];
    p = -grad;
    
    p0 = -grad_a0;
    p1 = -grad_a1;
    p2 = -grad_a2;
    
    while double(E(a(1) + lr * p0, a(2) + lr * p1, a(3) + lr * p2))...
            > double(E(a(1), a(2), a(3)) + lr * beta * grad' * p)
        
        lr = tau * lr;
        
    end
    
    % update ak
    a(1) = a(1) + lr * p0;
    a(2) = a(2) + lr * p1;
    a(3) = a(3) + lr * p2;
    
    step = abs(double(E(a(1), a(2), a(3))) - E_pre);
    
    
    E_record = [E_record; E_pre];
    distance_record = [distance_record; (a - optimal)' * (a - optimal)];
    
    count = count + 1;
    disp(a)
    
end

figure(1)
plot(10 * (1 : length(E_record)), real(E_record));
title('E vs Iteration');
xlabel('number of iteration');
ylabel('E value');

figure(2)
plot(10 * (1 : length(distance_record)), real(distance_record));
title('distance to optimal solution vs Iteration');
xlabel('number of iteration');
ylabel('distance');