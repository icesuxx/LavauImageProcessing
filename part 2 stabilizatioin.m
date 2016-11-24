clear all
close all
clc

Matrix = imread('54.png'); %first image

GrayI = rgb2gray(Matrix); %ทำภาพให้เป็นขาว-ดำ

[MatrixRow ,MatrixCol] = size(GrayI); %ขนาดของเมริกภาพ 100*1920

DoubleGrayI = mat2gray(GrayI); %แปลงเมตริกภาพขาวดำ จาก uint8 เปน double

%--------------%

%หา Gradient magnitude หรือความแรงของขอบภาพ

[Gx, Gy] = imgradientxy(GrayI); %Directional gradients along x-axis (horizontal) and y-axis (vertical)

[Gmag, Gdir] = imgradient(Gx, Gy); %returns the gradient magnitude and direction using directional gradients.

%--------------%

edge = 50; %ค่าขอบภาพกำหนดเอง

DivideRow = 245; %แบ่งภาพตาม row
DivideCol = 455; %แบ่งภาพตาม column

InsideMatrixRow = MatrixRow-(edge*2); %1080-(50*2)=980
InsideMatrixCol = MatrixCol-(edge*2); %1920-(50*2)=1820

SubmatrixRowSize = InsideMatrixRow/DivideRow; %980/245=4
SubmatrixColSize = InsideMatrixCol/DivideCol; %1820/455=4

a = 1;

for i = edge+1:DivideRow:MatrixRow-(edge+1) %เริ่มตั้งแต่ 51:245:1029
    for j = edge+1:DivideCol:MatrixCol-(edge+1) %51:455:1869
        
        Submatrix = Gmag((i:i+(DivideRow-1)),(j:j+DivideCol-1)); %แบ่งรูปจากเมตริก Gmag
        output = ['Submatrix' int2str(a) '=Submatrix'];
        eval(output);
        SubmatrixSum(a) = sum(sum(eval(['Submatrix' int2str(a)]))); %เก็บค่าของแต่ละเมตริกไว้
        MaxSubmatrixSum = max(SubmatrixSum);
        FindMaxSubmatrix = find(SubmatrixSum==max(SubmatrixSum)); %(เมตริก 10) ทำการหาค่า Gmag ที่มากที่สุด
        
        a = a+1;
        
    end
end

MaxSubmatrix1 = eval(['Submatrix' int2str(FindMaxSubmatrix)]); %เมตริกที่มีค่าGmagรวมทั้งเมตริกซ์แล้วมากที่สุด ตอบ เมตริก 10

MaxSubmatrix2 = Submatrix9; %เมตริกที่มีค่าGmagรวมทั้งเมตริกซ์แล้วมากอันดับที่ 2

%show image MaxSubmatrix1 and MaxSubmatrix2 in gray scale
figure,imagesc(MaxSubmatrix1);colormap(gray);axis off,axis image;
figure,imagesc(MaxSubmatrix2);colormap(gray);axis off,axis image;

%----------%

%size of MaxSubmatrix
[SubmatrixRow, SubmatrixCol] = size(MaxSubmatrix1); 

found = 0;
y = 0;

%----------%
%ทำการแบ่งภาพแบบเก็บพิกัด top-leftไว้ คือตัวแปร X
pp = 1;
for m = edge+1:DivideRow:(edge+1)+(DivideRow*floor(SubmatrixRowSize)-1)
    % for i = 3:2:5
    for n = edge+1:DivideCol:(edge+1)+(DivideCol*floor(SubmatrixColSize)-1)
        X(pp,1:2) = [m n];
        pp = pp + 1;
    end
end

%พิกัด top-left ของ MaxSubmatrix1 (เมตริกซ์ที่มีค่า Gmag มากที่สุด)
cor1 = [X(FindMaxSubmatrix,1),X(FindMaxSubmatrix,2)];
roww1 = cor1(1);
coll1 = cor1(2);

%พิกัด top-left ของ MaxSubmatrix2 (เมตริกซ์ที่มีค่า Gmag มากทอันดับ2)
cor2 = [X(9,1),X(9,2)];
roww2 = cor2(1);
coll2 = cor2(2);

fprintf('Top-Left coordinate Submatrix in Matrix(%d,%d)\n',roww1,coll1);
fprintf('Top-Left coordinate SubmatrixNew in Matrix(%d,%d)\n',roww2,coll2);
%----------%

%import รูปถัดไปเพื่อทำการหาค่าความ Stabilization
NewSnap = imread('55.png');

GrayI2 = rgb2gray(NewSnap); %ทำภาพให้เป็นขาว-ดำ

[Gx2, Gy2] = imgradientxy(GrayI2);

[Gmag2, Gdir2] = imgradient(Gx2, Gy2);

NewMatrix = Gmag2;

%----------%

[MatrixRowNew ,MatrixColNew] = size(NewMatrix);

foundNew = 0;
q1 = roww1; p1 = coll1;
q2 = roww2; p2 = coll2;
range = 40; %range ที่ใช้ในการแสกนหาเมตริกซ์ย่อยของรูปแรกในรูปที่2 คือ บน-ล่าง ซ้าย-ขวา
b = 1; mint = 1e15; cx = 0; dx = 0;

for c = -range : range
    for d = -range : range
        SubmatrixNewFirst = NewMatrix(q1+c:q1+c+(SubmatrixRow-1),p1+d:p1+d+(SubmatrixCol-1));
        SubmatrixNewSecond = NewMatrix(q2+c:q2+c+(SubmatrixRow-1),p2+d:p2+d+(SubmatrixCol-1));
        MatrixDiffFirst = abs(MaxSubmatrix1-SubmatrixNewFirst).^2;
        MatrixDiffSecond = abs(MaxSubmatrix2-SubmatrixNewSecond).^2;
        SumMatrixDiff = MatrixDiffFirst+MatrixDiffSecond;
        SumMatrixDiffM(b) = sum(sum(SumMatrixDiff));
        
        if SumMatrixDiffM(b) < mint
            mint = SumMatrixDiffM(b);
            cx = c+q1;
            dx = d+p1;
        end
        b = b+1;
    end
    % Move to next row
    q1 = q1+1;
    q2 = q2+1;
    % Reset column to initial
    p1 = coll1;
    p2 = coll2;
end

%หาเมตริกซ์ที่คำนวณผลต่างออกมากแล้วน้อยที่สุด
FindMinMatrixDiff = min(SumMatrixDiff);

FindMinLocaMatrixDiff = find(min(SumMatrixDiffM)==SumMatrixDiffM);

%บอกพิกัด top-left ของเมตริกซ์ผลตต่างน้อยที่สุดในรูปที่2
fprintf('Top-Left coordinate 2Submatrix in NewMatrix(%d,%d)\n',cx,dx);

%-------------%

%คำนวณออกมาว่าพิกัด top-left ของรูป2 ต่างกับรูปแรกเท่าไหร่  เอาไว้ดึงรูปกลับ หรือทำ stabilization
RowShift = cx - roww1;
ColShift = dx - coll1;

fprintf('Row Shift %d\n',RowShift);
fprintf('Column Shift %d\n',ColShift);
