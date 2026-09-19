clear;
clc;
close all;

%%Change the numLevels for different wavefunctions

%%parameters
numLevels = 7; 
hbar = 1;
m  = 1;
L = 1;

N = 1000;

dx = L/(N+1);

x = (1:N)' * dx; %%Divide into subintervals and transposes x

%%Kinetic energy operator
mainDiagonal = -2 * ones(N,1);
offDiagonal = 1 * ones(N-1,1);

D2 = (diag(mainDiagonal)...
    + diag(offDiagonal,1)...
    + diag(offDiagonal,-1)) / dx^2;

T = -(hbar^2/(2*m)) * D2;

%%Potential energy

V = zeros(N,1);

Vmatrix = diag(V);

%%Hamiltonian

H = T + Vmatrix;

%%Solving the eigenvalue problem

[psi,E] = eig(H);

E = diag(E);

%% Normalize wavefunctions
for n = 1:numLevels
    psi(:,n) = psi(:,n) / ...
        sqrt(trapz(x,abs(psi(:,n)).^2));
end

%% Analytical energies
n = (1:numLevels)';

E_exact = n.^2 * pi^2 / 2;

%% Display results
fprintf('n   Numerical E   Analytical E   Error\n');

for i = 1:numLevels
    error = abs(E(i)-E_exact(i));

    fprintf('%d   %.8f   %.8f   %.8e\n', ...
        i,E(i),E_exact(i),error);
end

%% Plot wavefunctions
figure;
hold on;

for n = 1:numLevels
    plot(x, psi(:,n) + E(n), ...
        'LineWidth', 1.5);
end

xlabel('x');
ylabel('\psi_n(x) + E_n');
title('Infinite Square Well: Numerical Wavefunctions');
grid on;