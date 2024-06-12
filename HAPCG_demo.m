% ***************************************************************************************************************************************
%  ***** HAPCG Alogrithm *******
%        [1]Y. Yao, Y. Zhang, Y. Wan, X. Liu, H. Guo. Heterologous Images Matching Considering Anisotropic Weighted Moment and 
%           Absolute Phase Orientation[J]. Geomatics and Information Science of Wuhan University, 2021, 46(11): 1727-1736. 
%           https://doi.org/10.13203/j.whugis20200702.
%        This is a simplified Code demo of the HAPCG algorithm.
%        Download website address of code and Images dataset:    https://skyearth.org/research
%        Public: Created by Yongxiang Yao in 2021/03/29.
%  ***************************************************************************************************************************************
clc
clear ;
close all;
warning('off');
%% 1 Import and display reference and image to be registered
file_image1= '.\Images\1DSMsets\pair1-2.jpg';
file_image2= '.\Images\1DSMsets\pair1-1.jpg';


image_1=imread(file_image1);
image_2=imread(file_image2);
[m,n,c1] = size(image_1);
[i,j,c2] = size(image_2);
if m ~= n || i~=j  || m ~= i 
    kk=max([m,n,i,j]);
    image_1 = imresize(image_1,[kk, kk]);
    image_2 = imresize(image_2,[kk, kk]);

end
save_path  = '.\outputs\1DSMsets\pair1r.jpg';
%% -----------Additive noise--------------- % SNR=0;
% SNR=-1;
% [image_1,noise]=Additive_noise(image_1,SNR);
% save_path  = '.\outputs\1DSMsets\add_pair1+0.jpg';

%% ------------Stripe noise---------------- % sigma=0.10;
% sigma=0.12;
% image_1 = Stripe_noise(image_1,sigma);
% save_path  = '.\outputs\1DSMsets\stripe_pair1+0p10.jpg';

%% 2  Setting of initial parameters 
% Key parameters:
K_weight=3;                       % 各向异性力矩图的加权值，处于（1~10），默认设置：3
Max=3;                            % Number of levels in scale space，默认设置：3
threshold = 0.4;                  % 特征点提取阈值，对SAR影像/强度图配色时，设置为：0.3；一般默认设置为：0.4
scale_value=2;                    % 尺度缩放比例值，默认设置：1.6
Path_Block=42;                    % 描述子邻域窗口大小， 默认设置：42；当需要更多特征点时，可以调大窗口。

%% 3 Anisotropic scale space
t1=clock;
disp('Start HAPCG algorithm processing, please waiting...');
tic;
[nonelinear_space_1]=HAPCG_nonelinear_space(image_1,Max,scale_value);      %-----------noise-free image
% [nonelinear_space_1]=HAPCG_nonelinear_space(image_1N,Max,scale_value);   %-----------noise image
[nonelinear_space_2]=HAPCG_nonelinear_space(image_2,Max,scale_value);
disp(['the cost time of constructinge anisotropic scale space：',num2str(toc),'秒']);

%% 4 constructing the anisotropic weighted moment map and  phase consistency orientation gradient
tic;
[harris_function_1,gradient_1,angle_1]=HAPCG_Gradient_Feature(nonelinear_space_1,Max,K_weight);
[harris_function_2,gradient_2,angle_2]=HAPCG_Gradient_Feature(nonelinear_space_2,Max,K_weight);
disp(['constructing the absolute phase consistency gradient map:',num2str(toc),'S']);

%% 5  feature point extraction
tic;
position_1=Harris_extreme(harris_function_1,gradient_1,angle_1,Max,threshold);
position_2=Harris_extreme(harris_function_2,gradient_2,angle_2,Max,threshold);
disp(['Feature point extraction time:  ',num2str(toc),' S']);

%% 6 Lop-Polar Descriptor Constrained by HAPCG
tic;
descriptors_1=HAPCG_Logpolar_descriptors(gradient_1,angle_1,position_1,Path_Block);                                     
descriptors_2=HAPCG_Logpolar_descriptors(gradient_2,angle_2,position_2,Path_Block); 
disp(['HAPCG Descriptor cost time:  ',num2str(toc),'S']); 

%% 7 Nearest matching    
disp('Nearest matching')
[indexPairs,~] = matchFeatures(descriptors_1.des,descriptors_2.des,'MaxRatio',1,'MatchThreshold', 10);
matchedPoints_1 = descriptors_1.locs(indexPairs(:, 1), :);
matchedPoints_2 = descriptors_2.locs(indexPairs(:, 2), :);
%% Outlier removal  
disp('Outlier removal')
[H,rmse]=FSC(matchedPoints_1,matchedPoints_2,'affine',3);
Y_=H*[matchedPoints_1(:,[1,2])';ones(1,size(matchedPoints_1,1))];
Y_(1,:)=Y_(1,:)./Y_(3,:);
Y_(2,:)=Y_(2,:)./Y_(3,:);
E=sqrt(sum((Y_(1:2,:)-matchedPoints_2(:,[1,2])').^2));
inliersIndex=E < 3;
clearedPoints1 = matchedPoints_1(inliersIndex, :);
clearedPoints2 = matchedPoints_2(inliersIndex, :);
uni1=[clearedPoints1(:,[1,2]),clearedPoints2(:,[1,2])];
[~,i,~]=unique(uni1,'rows','first');
inliersPoints1=clearedPoints1(sort(i)',:);
inliersPoints2=clearedPoints2(sort(i)',:);
[inliersPoints_1,inliersPoints_2] = BackProjection(inliersPoints1,inliersPoints2,scale_value);  % ---Projection to the original scale

disp('keypoints numbers of outlier removal: '); disp(size(inliersPoints_1,1));
disp(['RMSE of Matching results: ',num2str(rmse),'  pixel']);
figure; 
Matchresult=showMatchedFeatures(image_1, image_2, inliersPoints_1, inliersPoints_2, 'montage');  %-----------noise image
t2=clock;
disp(['Total time cost on HAPCG algorithm matching   :',num2str(etime(t2,t1)),' S']); 
correct_num=size(inliersPoints_1,1);

text(50,10,'HAPCG:','horiz','center','color','r','FontSize',15)
text(60,30,['Correct num:',num2str(correct_num)],'horiz','center','color','r','FontSize',15);
text(60,50,['RMSE:',num2str(rmse)],'horiz','center','color','r','FontSize',15);
text(60,70,['RT:',num2str(etime(t2,t1)),'s'],'horiz','center','color','r','FontSize',15);

set(gcf,'color',[1,1,1]) 
f=getframe(gcf);
imwrite(f.cdata,save_path)

    
