close all; clear all; clc;

syms i1 i2 i3 i4 i5 i6

%resistor values (ohms):
r = [10, 20, 2, 5, 25, 5];

%Circuit 1: 6 equations, 6 unknowns
%KCL equations:
% 0 = i1 - i2 - i6
% 0 = i2 - i3
% 0 = i3 - i4
% 0 = i4 - i5 + i6
%KVL equations:
% -200 = 10*i1 + 25*i5 + 5*i6
% 0 = 20*i2 + 2*i3 + 5*i4 - 5*i6

%Circuit 1 in a diagonally dominant matrix:
A1 = [1, -1, 0, 0, 0, -1; ...
      0, 20, 2, 5, 0, -5; ...
      0, 1, -1, 0, 0, 0; ...
      0, 0, 1, -1, 0, 0; ...
      10, 0, 0, 0, 25, 5; ...
      0, 0, 0, 1, -1, 1];
B1 = [0, 0, 0, 0, -200, 0];
%Circuit 1 equations (in order of matrix)
eq1 = i1 - i2 - i6 == 0;
eq2 = 20*i2 + 2*i3 + 5*i4 - 5*i6 == 0;
eq3 = i2 - i3 == 0;
eq4 = i3 - i4 == 0;
eq5 = 10*i1 + 25*i5 + 5*i6 == -200;
eq6 = i4 - i5 + i6 == 0;   

%Solve with Gauss Seidel - initial guesses = zeroes -- zeros(size(B))
convergence = 0.000001;
noRelaxationGS = GaussSeidel(A1, B1, zeros(size(B1)), 1.00, convergence);
relaxedGS = GaussSeidel(A1, B1, zeros(size(B1)), 1.02, convergence);
underRelaxedGS = GaussSeidel(A1, B1, zeros(size(B1)), 0.98, convergence);

%Compare results to matlab equation solver (m)
m = solve([eq1, eq2, eq3, eq4, eq5, eq6], [i1, i2, i3, i4, i5, i6]);
matlabSolution = double([m.i1, m.i2, m.i3, m.i4, m.i5, m.i6]);
%Solve using inverse
inverseSolution = (inv(A1)*B1.')';

%Display results - currents
fprintf("Circuit 1:\n")
colNames = {'i1', 'i2', 'i3', 'i4', 'i5', 'i6'};
I = [noRelaxationGS; relaxedGS; underRelaxedGS; matlabSolution; inverseSolution];
rowNames = {'No relaxation', 'Over relaxation', 'Under relaxation', 'Matlab solution', 'Inverse'};
circuit1I = table(I(:,1), I(:,2), I(:,3), I(:,4), I(:,5), I(:,6), 'VariableNames', colNames, 'RowNames', rowNames);
disp(circuit1I)

%Display results - voltages -> V=IR
colNames = {'v1', 'v2', 'v3', 'v4', 'v5', 'v6'};
V = [r.*noRelaxationGS; r.*relaxedGS; r.*underRelaxedGS; r.*matlabSolution; r.*inverseSolution];
circuit1V = table(V(:,1), V(:,2), V(:,3), V(:,4), V(:,5), V(:,6), 'VariableNames', colNames, 'RowNames', rowNames);
disp(circuit1V)