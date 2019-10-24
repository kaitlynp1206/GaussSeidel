close all; clear all; clc;

syms i1 i2 i3 i4 i5 

%resistor values (ohms):
r = [5, 10, 20, 15, 25];

%Circuit 3: 5 equations, 5 unknowns
%KCL equations:
% 0 = i1 - i2 - i4
% 0 = i2 + i3 - i5
%KVL equations:
% -80 = 5*i1 + 15*i4
% 0 = 10*i2 - 15*i4 + 25*i5
% -50 = 20*i3 + 25*i5

%Circuit 3 in a diagonally dominant matrix:
A3 = [1, -1, 0, -1, 0; ...
      0, 1, 1, 0, -1; ...
      0, 0, 20, 0, 25; ...
      5, 0, 0, 15, 0; ...
      0, 10, 0, -15, 25];
B3 = [0, 0, -50, -80, 0];
%Circuit 3 equations (in order of matrix)
eq1 = i1 - i2 - i4 == 0;
eq2 = i2 + i3 - i5 == 0;
eq3 = 20*i3 + 25*i5 == -50;
eq4 = 5*i1 + 15*i4 == -80;
eq5 = 10*i2 - 15*i4 + 25*i5 == 0;

%Solve with Gauss Seidel - initial guesses = zeroes -- zeros(size(B))
convergence = 0.000001;
noRelaxationGS = GaussSeidel(A3, B3, zeros(size(B3)), 1.00, convergence);
relaxedGS = GaussSeidel(A3, B3, zeros(size(B3)), 1.02, convergence);
underRelaxedGS = GaussSeidel(A3, B3, zeros(size(B3)), 0.98, convergence);

%Compare results to matlab equation solver (m)
m = solve([eq1, eq2, eq3, eq4, eq5], [i1, i2, i3, i4, i5]);
matlabSolution = double([m.i1, m.i2, m.i3, m.i4, m.i5]);
%Solve using inverse
inverseSolution = (inv(A3)*B3.')';

%Display results - currents
fprintf("Circuit 3:\n")
colNames = {'i1', 'i2', 'i3', 'i4', 'i5'};
I = [noRelaxationGS; relaxedGS; underRelaxedGS; matlabSolution; inverseSolution];
rowNames = {'No relaxation', 'Over relaxation', 'Under relaxation', 'Matlab solution', 'Inverse'};
circuit3I = table(I(:,1), I(:,2), I(:,3), I(:,4), I(:,5), 'VariableNames', colNames, 'RowNames', rowNames);
disp(circuit3I)

%Display results - voltages -> V=IR
fprintf("Circuit 3:\n")
colNames = {'v1', 'v2', 'v3', 'v4', 'v5'};
V = [r.*noRelaxationGS; r.*relaxedGS; r.*underRelaxedGS; r.*matlabSolution; r.*inverseSolution];
circuit3V = table(V(:,1), V(:,2), V(:,3), V(:,4), V(:,5), 'VariableNames', colNames, 'RowNames', rowNames);
disp(circuit3V)