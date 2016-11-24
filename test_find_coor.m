clc
clear all
close all

Matrix = [10 11 12 13 14 15 16 17 18 19;
    20 21 22 23 24 25 26 27 28 29;
    30 31 32 33 34 35 36 37 38 39;
    40 41 42 43 44 45 46 47 48 49;
    50 51 52 53 54 55 56 57 58 59;
    60 61 62 63 64 65 66 67 68 69;
    70 71 72 73 74 75 76 77 78 79;
    80 81 82 83 84 85 86 87 88 89;
    90 91 92 93 94 95 96 97 98 99];

[MatrixRow ,MatrixCol] = size(Matrix);
% MatrixRow = 1080;
% MatrixCol = 1920;

% edge = 50;
% row = 245;
% col = 455;

edge = 2;
row = 2;
col = 2;

InsideMatrixRow = MatrixRow-(edge*2);
InsideMatrixCol = MatrixCol-(edge*2);

SubMatrixRowSize = InsideMatrixRow/row;
SubMatrixColSize = InsideMatrixCol/col;

pp = 1;
for i=edge+1:row:(edge+1)+(row*floor(SubMatrixRowSize)-1)
% for i=3:2:5
for j=edge+1:col:(edge+1)+(col*floor(SubMatrixColSize)-1)
%     for j=3:2:7
        X(pp,1:2) = [i j];
        pp = pp + 1;
    end
end

% cor = [X(16,1),X(16,2)];
cor = [X(6,1),X(6,2)];

roww = cor(1);
coll = cor(2);
