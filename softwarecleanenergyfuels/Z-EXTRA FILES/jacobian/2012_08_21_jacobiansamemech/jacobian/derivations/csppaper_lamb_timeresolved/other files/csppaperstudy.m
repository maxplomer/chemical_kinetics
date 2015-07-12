clc;clear;
syms a b c d e f g h i;
A = [a b;c d]
% A = [1 2;
%      3 4]
B = inv(A)

B(1,1)*A(1,1)+B(1,2)*A(2,1)
B(1,1)*A(1,2)+B(1,2)*A(2,2)
