clear all
close all
clc

%% frame rate 10 เอาไว้หา median เพราะถ้าเอา famerate ต่ำๆมาหาจะเกิด noise เยอะ
I1 = imread('53.png');
I2 = imread('63.png');
I3 = imread('73.png');
I4 = imread('83.png');
I5 = imread('93.png');
I6 = imread('103.png');
I7 = imread('113.png');

%% frame rate 1 เอาไว้หาจำนวนลูกน้ำ และประมวลผลภาพ
I8 = imread('54.png');
I9 = imread('55.png');
I10 = imread('56.png');
I11 = imread('57.png');
I12 = imread('58.png');
I13 = imread('59.png');
I14 = imread('60.png');

%%ทำภาพที่ import มาให้เป็น gray scale
GrayI1 = rgb2gray(I1);
GrayI2 = rgb2gray(I2);
GrayI3 = rgb2gray(I3);
GrayI4 = rgb2gray(I4);
GrayI5 = rgb2gray(I5);
GrayI6 = rgb2gray(I6);
GrayI7 = rgb2gray(I7);

GrayI8 = rgb2gray(I8);
GrayI9 = rgb2gray(I9);
GrayI10 = rgb2gray(I10);
GrayI11 = rgb2gray(I11);
GrayI12 = rgb2gray(I12);
GrayI13 = rgb2gray(I13);
GrayI14 = rgb2gray(I14);

edge = 50; %ค่าขอบภาพ

a = 1;

[MatrixRow, MatrixCol] = size(GrayI1);

%% trip image (cut edge) ทำทั้งหมด 14 รูปเลย

for i = 1:14
    
    DoubleGrayImage = mat2gray(eval(['GrayI' int2str(i)]));
    
    output = ['DoubleGrayI' int2str(i) '=DoubleGrayImage'];
    eval(output);
    %เอาแค่ตัวรูปข้างในถัดจาก edge เข้าไป
    InsideMatrix = DoubleGrayImage(edge+1:MatrixRow-(edge+1),...
        edge+1:MatrixCol-(edge+1));  %แบ่งรูปจากเมตริก Gmag
    output2 = ['InsideMatrixI' int2str(i) '=InsideMatrix'];
    eval(output2);
    
    GrayInsideMatrix = mat2gray(eval(['InsideMatrixI' int2str(i)]));
    output3 = ['GrayInsideMatrix' int2str(i) '=GrayInsideMatrix'];
    eval(output3);
    
end

%% 1.find Background (Median)
% 2.divide image from Background
% 3.find FirstThreshold (Mean and Max 2.)
% 4.vary = 0.75

%เอาเมตริกซ์มาวางซ้อนๆกันไว้
z = cat(3,InsideMatrixI1...
    ,InsideMatrixI2,InsideMatrixI3,InsideMatrixI4,InsideMatrixI5,InsideMatrixI6,InsideMatrixI7);

DoubleZ = mat2gray(z);

%หาค่า median ของรูป 10 framrate ทุกรูปที่วางซ้อนทับไว้ เรียกว่ารูป background
BackgroundIm = median(DoubleZ,3);

for j = 1:14
    
	%เอารุป 1 framrateทุกๆรูปมาลบกับรูป background แล้วจะได้รูปใหม่ที่มีแต่ลูกน้ำ หรือสิ่งที่เคลื่อนไหว
    BgDivideInsideMatrix =  abs(eval(['GrayInsideMatrix' int2str(j)])-BackgroundIm);
    output4 = ['BgDivideInsideMatrix' int2str(j) '=BgDivideInsideMatrix'];
    eval(output4);
    
	%หาค่า mean ของรูปใหม่หลังบลกับ background
    MeanInsideMatrix = mean2(eval(['BgDivideInsideMatrix' int2str(j)]));
    output5 = ['MeanInsideMatrix' int2str(j) '=MeanInsideMatrix'];
    eval(output5);
    
	%หาค่า max ของ pixel จากรูปที่ลบ background แล้ว หรือค่าสีที่สว่างที่สุด
    [MaxInsideMatrix, Location] = max(eval(['BgDivideInsideMatrix' int2str(j)]));
    output5 = ['MaxInsideMatrix' int2str(j) '=MaxInsideMatrix'];
    eval(output5);
    
	%กำหนดค่านี้เองว่าให้แวรี่ = 0.7
    vary = 0.7;
    
	%FirstThreshold เกิดจากสูตร  ((vary*MeanInsideMatrix)+(1-vary)*MeanInsideMatrix))
	FirstThreshold = (vary*(eval(['MeanInsideMatrix' int2str(j)])))+((1-vary)*(eval(['MaxInsideMatrix' int2str(j)])));
    output6 = ['FirstThreshold' int2str(j) '=FirstThreshold'];
    eval(output6);
    
end

%% Transform Image after divide from Background to 2 value
% 1. 255 when > FirstThreshold
% 2. 0 when <,= FirstThreshold

[RowInsideMatrix, ColInsideMatrix] = size(BgDivideInsideMatrix1);

%สร้างเมตริก 1 ที่มีขนาดเท่ากับรูปหลังลบ bachground ไว้ใส่ข้อมูล
MatrixTrans = ones(RowInsideMatrix,ColInsideMatrix);

for k = 1:14
    aaa = eval(['BgDivideInsideMatrix' int2str(k)]);
    for c = 1:numel(aaa)
        if aaa(c) > eval(['FirstThreshold' int2str(k)])%ถ้าค่าใน pixel นั้นมากกว่า FirstThreshold ให้มีค่าเปน 255 หรือสีขาว
            MatrixTrans(c) = 255;
        else  %ถ้าค่าใน pixel นั้นน้อยกว่า FirstThreshold ให้มีค่าเปน 0 หรือสีดำ
            MatrixTrans(c) = 0;
        end
    end
    
    output7 = ['MatrixTrans' int2str(k) '=MatrixTrans'];
    eval(output7);
    
end

for l = 1:14
    
    BlackWhiteMatrix = im2bw(eval(['MatrixTrans' int2str(l)]),1);
    output8 = ['BlackWhiteMatrix' int2str(l) '=BlackWhiteMatrix'];
    eval(output8);
    
end

% figure, imagesc(BlackWhiteMatrix1);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix2);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix3);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix4);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix5);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix6);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix7);colormap(gray);axis off,axis image;

figure, imagesc(BackgroundIm);colormap(gray);axis off,axis image;

% figure, imagesc(BlackWhiteMatrix8);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix9);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix10);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix11);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix12);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix13);colormap(gray);axis off,axis image;
% figure, imagesc(BlackWhiteMatrix14);colormap(gray);axis off,axis image;

%% find BlackWhiteMatrixDetail :

for n = 8:14 %เอาเฉพาะรูปที่ framerate 1 มาหาจำนนลูกน้ำ
    
    BlackWhiteMatrixDetail = bwconncomp(eval(['BlackWhiteMatrix' int2str(n)])); %bwconncomp เธ?เธทเธญเธ?เธณเธชเธฑเน?เธ?เน€เธ?เธ?เธฒเธฐเธ?เธญเธ? MATLAB connectivity,ImageSize,NumObjects,PixelIdxList
    output9 = ['BlackWhiteMatrixDetail' int2str(n) '=BlackWhiteMatrixDetail'];
    eval(output9);
    
end

%บอกจำนวน object ในแต่ละรูป
NumOfObject8 = BlackWhiteMatrixDetail8.NumObjects;
NumOfObject9 = BlackWhiteMatrixDetail9.NumObjects;
NumOfObject10 = BlackWhiteMatrixDetail10.NumObjects;
NumOfObject11 = BlackWhiteMatrixDetail11.NumObjects;
NumOfObject12 = BlackWhiteMatrixDetail12.NumObjects;
NumOfObject13 = BlackWhiteMatrixDetail13.NumObjects;
NumOfObject14 = BlackWhiteMatrixDetail14.NumObjects;

%% put square frame over object

for m = 8:14
    
    BW = im2bw(eval(['BlackWhiteMatrix' int2str(m)]));
    output10 = ['BW' int2str(m) '=BW'];
    eval(output10);
    
    StatsBW = regionprops(eval(['BW' int2str(m)])); %regionprops คือคำสั่งเฉพาะของ MATALB บอกค่า Area, Centroid, BoundingBox
    output11 = ['StatsBW' int2str(m) '=StatsBW'];
    eval(output11);
    
    StatsNotBW = regionprops(eval(['BW' int2str(m)]));
    output12 = ['StatsNotBW' int2str(m) '=StatsNotBW'];
    eval(output12);
end
%% ภาพที่ 1 : 54
for l8 = 1:NumOfObject8 %วนไปทุก object ในรูปที่54
    
    [WhiteObjectSize8, nan] = size(cell2mat(BlackWhiteMatrixDetail8.PixelIdxList(l8))); %เมตริกเก็บขนาดของแต่ละ object
    
    WhiteObject8(l8) = WhiteObjectSize8;
    
    SortWhiteObject8 = sort(WhiteObject8); %เรียงจากขนาดน้อยไปมาก///
    
    MeanSizeOfObject8 = mean(WhiteObject8); %หาค่า mean ของขนาด object///
    
    SdSizeOfObject8 = std(WhiteObject8); % หาค่า standard deviation ของขนาด object///
    
    [Nan1, NumWhiteObject8] = size(find(WhiteObject8));
    
    MinimunSizeMatrix8 = find(WhiteObject8 < 22); %pixel สีขาวติดกันอย่างน้อย 22 pixel
    MaximunSizeMatrix8 = find(WhiteObject8 > 134); %pixel สีขาวติดกันอย่างมาก 134 pixel
    
    FindWhiteObject8 = find(WhiteObject8>22 & WhiteObject8<134); %กรองขนาดของ object
    
    [Nan2, HalfMeanSizeOfObject8]= size(MinimunSizeMatrix8); 
    [Nan3, MeanPlusSdSizeOfObject8]= size(MaximunSizeMatrix8);
    
    MinimumWhiteObject8 = NumWhiteObject8 - HalfMeanSizeOfObject8;
    MaximumWhiteObject8 = NumWhiteObject8 - MeanPlusSdSizeOfObject8;
    
    NumOfLuava8 = (MinimumWhiteObject8 + MaximumWhiteObject8) - NumWhiteObject8; %บอกจำนวนลูกนำในรูปนั้น
    
    l8 = l8+1;
end

WhiteObjectPassFilter8 = StatsBW8(FindWhiteObject8);
[statsRow1, statsCol1] = size(WhiteObjectPassFilter8);

for o1 = 1:statsRow1
    WhitePixilInSquare1(o1) = WhiteObjectPassFilter8(o1).Area;

    CentroidRowOfBlock1(o1) = WhiteObjectPassFilter8(o1).Centroid(1);
    CentroidColOfBlock1(o1) = WhiteObjectPassFilter8(o1).Centroid(2);

    x_widthOfBlock1(o1) = WhiteObjectPassFilter8(o1).BoundingBox(3);
    y_widthOfBlock1(o1) = WhiteObjectPassFilter8(o1).BoundingBox(4);

    Aspect1(o1) = min(x_widthOfBlock1(o1),y_widthOfBlock1(o1))/max(x_widthOfBlock1(o1),y_widthOfBlock1(o1));

    PixelInSquare1(o1) = x_widthOfBlock1(o1)*y_widthOfBlock1(o1);

    Occupancy1(o1) = WhitePixilInSquare1(o1)/PixelInSquare1(o1);

end
figure,imshow(BW8);
hold on;
for p1 = 1:numel(WhiteObjectPassFilter8)
    rec1 = rectangle('Position', WhiteObjectPassFilter8(p1).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen1 = WhiteObjectPassFilter8(p1).Centroid;
    text1 = text(cen1(1), cen1(2), sprintf('%d',p1),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class1 = xlsread('Class2.xlsx','Sheet2','C1:C61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy1,Aspect1,Class1,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');

%% เธ เธฒเธ?เธ—เธตเน? 2 : 55
for l9 = 1:NumOfObject9
    
    [WhiteObjectSize9, nan] = size(cell2mat(BlackWhiteMatrixDetail9.PixelIdxList(l9)));
    
    WhiteObject9(l9) = WhiteObjectSize9;
    
    SortWhiteObject9 = sort(WhiteObject9);
    
    MeanSizeOfObject9 = mean(WhiteObject9);
    
    SdSizeOfObject9 = std(WhiteObject9);
    
    [Nan1, NumWhiteObject9] = size(find(WhiteObject9));
    
    MinimunSizeMatrix9 = find(WhiteObject9 < 22);
    MaximunSizeMatrix9 = find(WhiteObject9 > 134);
    
    FindWhiteObject9 = find(WhiteObject9>22 & WhiteObject9<134);
    
    [Nan2, HalfMeanSizeOfObject9]= size(MinimunSizeMatrix9);
    [Nan3, MeanPlusSdSizeOfObject9]= size(MaximunSizeMatrix9);
    
    MinimumWhiteObject9 = NumWhiteObject9 - HalfMeanSizeOfObject9;
    MaximumWhiteObject9 = NumWhiteObject9 - MeanPlusSdSizeOfObject9;
    
    NumOfLuava9 = (MinimumWhiteObject9 + MaximumWhiteObject9) - NumWhiteObject9;
    
    l9 = l9+1;
end

WhiteObjectPassFilter9 = StatsBW9(FindWhiteObject9);
[statsRow2, statsCol2] = size(WhiteObjectPassFilter9);

for o2 = 1:statsRow2
    WhitePixilInSquare2(o2) = WhiteObjectPassFilter9(o2).Area;

    CentroidRowOfBlock2(o2) = WhiteObjectPassFilter9(o2).Centroid(1);
    CentroidColOfBlock2(o2) = WhiteObjectPassFilter9(o2).Centroid(2);

    x_widthOfBlock2(o2) = WhiteObjectPassFilter9(o2).BoundingBox(3);
    y_widthOfBlock2(o2) = WhiteObjectPassFilter9(o2).BoundingBox(4);

    Aspect2(o2) = min(x_widthOfBlock2(o2),y_widthOfBlock2(o2))/max(x_widthOfBlock2(o2),y_widthOfBlock2(o2));

    PixelInSquare2(o2) = x_widthOfBlock2(o2)*y_widthOfBlock2(o2);

    Occupancy2(o2) = WhitePixilInSquare2(o2)/PixelInSquare2(o2);

end
figure,imshow(BW9);
hold on;
for p2 = 1:numel(WhiteObjectPassFilter9)
    rec2 = rectangle('Position', WhiteObjectPassFilter9(p2).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen2 = WhiteObjectPassFilter9(p2).Centroid;
    text2 = text(cen2(1), cen2(2), sprintf('%d',p2),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class2 = xlsread('Class2.xlsx','Sheet2','F1:F61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy2,Aspect2,Class2,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');
%% เธ เธฒเธ?เธ—เธตเน? 3 : 56
for l10 = 1:NumOfObject10
    
    [WhiteObjectSize10, nan] = size(cell2mat(BlackWhiteMatrixDetail10.PixelIdxList(l10)));
    
    WhiteObject10(l10) = WhiteObjectSize10;
    
    SortWhiteObject10 = sort(WhiteObject10);
    
    MeanSizeOfObject10 = mean(WhiteObject10);
    
    SdSizeOfObject10 = std(WhiteObject10);
    
    [Nan1, NumWhiteObject10] = size(find(WhiteObject10));
    
    MinimunSizeMatrix10 = find(WhiteObject10 < 22);
    MaximunSizeMatrix10 = find(WhiteObject10 > 134);
    
    FindWhiteObject10 = find(WhiteObject10>22 & WhiteObject10<134);
    
    [Nan3, HalfMeanSizeOfObject10]= size(MinimunSizeMatrix10);
    [Nan3, MeanPlusSdSizeOfObject10]= size(MaximunSizeMatrix10);
    
    MinimumWhiteObject10 = NumWhiteObject10 - HalfMeanSizeOfObject10;
    MaximumWhiteObject10 = NumWhiteObject10 - MeanPlusSdSizeOfObject10;
    
    NumOfLuava10 = (MinimumWhiteObject10 + MaximumWhiteObject10) - NumWhiteObject10;
    
    l10 = l10+1;
end

WhiteObjectPassFilter10 = StatsBW10(FindWhiteObject10);
[statsRow3, statsCol3] = size(WhiteObjectPassFilter10);

for o3 = 1:statsRow3
    WhitePixilInSquare3(o3) = WhiteObjectPassFilter10(o3).Area;

    CentroidRowOfBlock3(o3) = WhiteObjectPassFilter10(o3).Centroid(1);
    CentroidColOfBlock3(o3) = WhiteObjectPassFilter10(o3).Centroid(2);

    x_widthOfBlock3(o3) = WhiteObjectPassFilter10(o3).BoundingBox(3);
    y_widthOfBlock3(o3) = WhiteObjectPassFilter10(o3).BoundingBox(4);

    Aspect3(o3) = min(x_widthOfBlock3(o3),y_widthOfBlock3(o3))/max(x_widthOfBlock3(o3),y_widthOfBlock3(o3));

    PixelInSquare3(o3) = x_widthOfBlock3(o3)*y_widthOfBlock3(o3);

    Occupancy3(o3) = WhitePixilInSquare3(o3)/PixelInSquare3(o3);

end
figure,imshow(BW10);
hold on;
for p3 = 1:numel(WhiteObjectPassFilter10)
    rec3 = rectangle('Position', WhiteObjectPassFilter10(p3).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen3 = WhiteObjectPassFilter10(p3).Centroid;
    text3 = text(cen3(1), cen3(2), sprintf('%d',p3),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class3 = xlsread('Class2.xlsx','Sheet2','I1:I61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy3,Aspect3,Class3,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');
%% เธ เธฒเธ?เธ—เธตเน? 4 : 57
for l11 = 1:NumOfObject11
    
    [WhiteObjectSize11, nan] = size(cell2mat(BlackWhiteMatrixDetail11.PixelIdxList(l11)));
    
    WhiteObject11(l11) = WhiteObjectSize11;
    
    SortWhiteObject11 = sort(WhiteObject11);
    
    MeanSizeOfObject11 = mean(WhiteObject11);
    
    SdSizeOfObject11 = std(WhiteObject11);
    
    [Nan1, NumWhiteObject11] = size(find(WhiteObject11));
    
    MinimunSizeMatrix11 = find(WhiteObject11 < 22);
    MaximunSizeMatrix11 = find(WhiteObject11 > 134);
    
    FindWhiteObject11 = find(WhiteObject11>22 & WhiteObject11<134);
    
    [Nan4, HalfMeanSizeOfObject11]= size(MinimunSizeMatrix11);
    [Nan4, MeanPlusSdSizeOfObject11]= size(MaximunSizeMatrix11);
    
    MinimumWhiteObject11 = NumWhiteObject11 - HalfMeanSizeOfObject11;
    MaximumWhiteObject11 = NumWhiteObject11 - MeanPlusSdSizeOfObject11;
    
    NumOfLuava11 = (MinimumWhiteObject11 + MaximumWhiteObject11) - NumWhiteObject11;
    
    l11 = l11+1;
end

WhiteObjectPassFilter11 = StatsBW11(FindWhiteObject11);
[statsRow4, statsCol4] = size(WhiteObjectPassFilter11);

for o4 = 1:statsRow4
    WhitePixilInSquare4(o4) = WhiteObjectPassFilter11(o4).Area;

    CentroidRowOfBlock4(o4) = WhiteObjectPassFilter11(o4).Centroid(1);
    CentroidColOfBlock4(o4) = WhiteObjectPassFilter11(o4).Centroid(2);

    x_widthOfBlock4(o4) = WhiteObjectPassFilter11(o4).BoundingBox(3);
    y_widthOfBlock4(o4) = WhiteObjectPassFilter11(o4).BoundingBox(4);

    Aspect4(o4) = min(x_widthOfBlock4(o4),y_widthOfBlock4(o4))/max(x_widthOfBlock4(o4),y_widthOfBlock4(o4));

    PixelInSquare4(o4) = x_widthOfBlock4(o4)*y_widthOfBlock4(o4);

    Occupancy4(o4) = WhitePixilInSquare4(o4)/PixelInSquare4(o4);

end
figure,imshow(BW11);
hold on;
for p4 = 1:numel(WhiteObjectPassFilter11)
    rec4 = rectangle('Position', WhiteObjectPassFilter11(p4).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen4 = WhiteObjectPassFilter11(p4).Centroid;
    text4 = text(cen4(1), cen4(2), sprintf('%d',p4),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class4 = xlsread('Class2.xlsx','Sheet2','L1:L61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy4,Aspect4,Class4,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');
%% เธ เธฒเธ?เธ—เธตเน? 5 : 58
for l12 = 1:NumOfObject12
    
    [WhiteObjectSize12, nan] = size(cell2mat(BlackWhiteMatrixDetail12.PixelIdxList(l12)));
    
    WhiteObject12(l12) = WhiteObjectSize12;
    
    SortWhiteObject12 = sort(WhiteObject12);
    
    MeanSizeOfObject12 = mean(WhiteObject12);
    
    SdSizeOfObject12 = std(WhiteObject12);
    
    [Nan1, NumWhiteObject12] = size(find(WhiteObject12));
    
    MinimunSizeMatrix12 = find(WhiteObject12 < 22);
    MaximunSizeMatrix12 = find(WhiteObject12 > 134);
    
    FindWhiteObject12 = find(WhiteObject12>22 & WhiteObject12<134);
    
    [Nan4, HalfMeanSizeOfObject12]= size(MinimunSizeMatrix12);
    [Nan4, MeanPlusSdSizeOfObject12]= size(MaximunSizeMatrix12);
    
    MinimumWhiteObject12 = NumWhiteObject12 - HalfMeanSizeOfObject12;
    MaximumWhiteObject12 = NumWhiteObject12 - MeanPlusSdSizeOfObject12;
    
    NumOfLuava12 = (MinimumWhiteObject12 + MaximumWhiteObject12) - NumWhiteObject12;
    
    l12 = l12+1;
end

WhiteObjectPassFilter12 = StatsBW12(FindWhiteObject12);
[statsRow5, statsCol5] = size(WhiteObjectPassFilter12);

for o5 = 1:statsRow5
    WhitePixilInSquare5(o5) = WhiteObjectPassFilter12(o5).Area;

    CentroidRowOfBlock5(o5) = WhiteObjectPassFilter12(o5).Centroid(1);
    CentroidColOfBlock5(o5) = WhiteObjectPassFilter12(o5).Centroid(2);

    x_widthOfBlock5(o5) = WhiteObjectPassFilter12(o5).BoundingBox(3);
    y_widthOfBlock5(o5) = WhiteObjectPassFilter12(o5).BoundingBox(4);

    Aspect5(o5) = min(x_widthOfBlock5(o5),y_widthOfBlock5(o5))/max(x_widthOfBlock5(o5),y_widthOfBlock5(o5));

    PixelInSquare5(o5) = x_widthOfBlock5(o5)*y_widthOfBlock5(o5);

    Occupancy5(o5) = WhitePixilInSquare5(o5)/PixelInSquare5(o5);

end
figure,imshow(BW12);
hold on;
for p5 = 1:numel(WhiteObjectPassFilter12)
    rec5 = rectangle('Position', WhiteObjectPassFilter12(p5).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen5 = WhiteObjectPassFilter12(p5).Centroid;
    text5 = text(cen5(1), cen5(2), sprintf('%d',p5),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class5 = xlsread('Class2.xlsx','Sheet2','O1:O61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy5,Aspect5,Class5,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');
%% เธ เธฒเธ?เธ—เธตเน? 6 : 59
for l13 = 1:NumOfObject13
    
    [WhiteObjectSize13, nan] = size(cell2mat(BlackWhiteMatrixDetail13.PixelIdxList(l13)));
    
    WhiteObject13(l13) = WhiteObjectSize13;
    
    SortWhiteObject13 = sort(WhiteObject13);
    
    MeanSizeOfObject13 = mean(WhiteObject13);
    
    SdSizeOfObject13 = std(WhiteObject13);
    
    [Nan1, NumWhiteObject13] = size(find(WhiteObject13));
    
    MinimunSizeMatrix13 = find(WhiteObject13 < 22);
    MaximunSizeMatrix13 = find(WhiteObject13 > 134);
    
    FindWhiteObject13 = find(WhiteObject13>22 & WhiteObject13<134);
    
    [Nan4, HalfMeanSizeOfObject13]= size(MinimunSizeMatrix13);
    [Nan4, MeanPlusSdSizeOfObject13]= size(MaximunSizeMatrix13);
    
    MinimumWhiteObject13 = NumWhiteObject13 - HalfMeanSizeOfObject13;
    MaximumWhiteObject13 = NumWhiteObject13 - MeanPlusSdSizeOfObject13;
    
    NumOfLuava13 = (MinimumWhiteObject13 + MaximumWhiteObject13) - NumWhiteObject13;
    
    l13 = l13+1;
end

WhiteObjectPassFilter13 = StatsBW13(FindWhiteObject13);
[statsRow6, statsCol6] = size(WhiteObjectPassFilter13);

for o6 = 1:statsRow6
    WhitePixilInSquare6(o6) = WhiteObjectPassFilter13(o6).Area;

    CentroidRowOfBlock6(o6) = WhiteObjectPassFilter13(o6).Centroid(1);
    CentroidColOfBlock5(o6) = WhiteObjectPassFilter13(o6).Centroid(2);

    x_widthOfBlock6(o6) = WhiteObjectPassFilter13(o6).BoundingBox(3);
    y_widthOfBlock6(o6) = WhiteObjectPassFilter13(o6).BoundingBox(4);

    Aspect6(o6) = min(x_widthOfBlock6(o6),y_widthOfBlock6(o6))/max(x_widthOfBlock6(o6),y_widthOfBlock6(o6));

    PixelInSquare6(o6) = x_widthOfBlock6(o6)*y_widthOfBlock6(o6);

    Occupancy6(o6) = WhitePixilInSquare6(o6)/PixelInSquare6(o6);

end
figure,imshow(BW13);
hold on;
for p6 = 1:numel(WhiteObjectPassFilter13)
    rec6 = rectangle('Position', WhiteObjectPassFilter13(p6).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen6 = WhiteObjectPassFilter13(p6).Centroid;
    text6 = text(cen6(1), cen6(2), sprintf('%d',p6),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
Class6 = xlsread('Class2.xlsx','Sheet2','R1:R61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy6,Aspect6,Class6,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');

%% เธ เธฒเธ?เธ—เธตเน? 7 : 60
for l14 = 1:NumOfObject14
    
    [WhiteObjectSize14, nan] = size(cell2mat(BlackWhiteMatrixDetail14.PixelIdxList(l14)));
    
    WhiteObject14(l14) = WhiteObjectSize14;
    
    SortWhiteObject14 = sort(WhiteObject14);
    
    MeanSizeOfObject14 = mean(WhiteObject14);
    
    SdSizeOfObject14 = std(WhiteObject14);
    
    [Nan1, NumWhiteObject14] = size(find(WhiteObject14));
    
    MinimunSizeMatrix14 = find(WhiteObject14 < 22);
    MaximunSizeMatrix14 = find(WhiteObject14 > 134);
    
    FindWhiteObject14 = find(WhiteObject14>22 & WhiteObject14<134);
    
    [Nan3, HalfMeanSizeOfObject14]= size(MinimunSizeMatrix14);
    [Nan3, MeanPlusSdSizeOfObject14]= size(MaximunSizeMatrix14);
    
    MinimumWhiteObject14 = NumWhiteObject14 - HalfMeanSizeOfObject14;
    MaximumWhiteObject14 = NumWhiteObject14 - MeanPlusSdSizeOfObject14;
    
    NumOfLuava14 = (MinimumWhiteObject14 + MaximumWhiteObject14) - NumWhiteObject14;
    
    l14 = l14+1;
end

WhiteObjectPassFilter14 = StatsBW14(FindWhiteObject14);
[statsRow7, statsCol7] = size(WhiteObjectPassFilter14);

for o7 = 1:statsRow7
    WhitePixilInSquare7(o7) = WhiteObjectPassFilter14(o7).Area;

    CentroidRowOfBlock7(o7) = WhiteObjectPassFilter14(o7).Centroid(1);
    CentroidColOfBlock7(o7) = WhiteObjectPassFilter14(o7).Centroid(2);

    x_widthOfBlock7(o7) = WhiteObjectPassFilter14(o7).BoundingBox(3);
    y_widthOfBlock7(o7) = WhiteObjectPassFilter14(o7).BoundingBox(4);

    Aspect7(o7) = min(x_widthOfBlock7(o7),y_widthOfBlock7(o7))/max(x_widthOfBlock7(o7),y_widthOfBlock7(o7));

    PixelInSquare7(o7) = x_widthOfBlock7(o7)*y_widthOfBlock7(o7);

    Occupancy7(o7) = WhitePixilInSquare7(o7)/PixelInSquare7(o7);

end
figure,imshow(BW14);
hold on;
for p7 = 1:numel(WhiteObjectPassFilter14)
    rec7 = rectangle('Position', WhiteObjectPassFilter14(p7).BoundingBox, ...
        'Linewidth', 2, 'EdgeColor', 'r', 'LineStyle', '-');
    cen7 = WhiteObjectPassFilter14(p7).Centroid;
    text7 = text(cen7(1), cen7(2), sprintf('%d',p7),...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom',...
        'Color','red');
end
%%
Class7 = xlsread('Class2.xlsx','Sheet2','U1:U61');
hold on;
% plot(CentroidRowOfBlock7, CentroidColOfBlock7, 'g+', 'MarkerSize', 20, 'LineWidth', 1);
% hold off;

figure, gscatter(Occupancy7,Aspect7,Class7,'rgb','osd'),axis([0 1.5 0 1.5]);
xlabel('Occupancy');
ylabel('Aspect');

toc;
