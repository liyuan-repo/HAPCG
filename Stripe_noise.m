function x_stripe = Stripe_noise(img1,sigma)


%-------------------- 乘性噪声----------------
% img=rgb2gray(img1);
img1=double(img1);
xorig = img1/255;
% xorig=double(img1);
[m,n,c] = size(xorig);
% u = sigma^2 *randn(m,n,c);  %---------彩色条纹
u = sigma^2 *randn(m,n);
uu(:,:,1)=u;
uu(:,:,2)=u;
uu(:,:,3)=u;                  %---------黑白条纹
kk=zeros(m,n,c);
for i=1:c
    kk(:,:,i)=xorig(:,:,i)*uu(:,:,i);
end
% kk=bsxfun(@times,u,xorig); % 这是元素乘积
xin = xorig+kk;
% xin = xorig + u * xorig;   % u与xorig之间注意区分元素相乘和矩阵相乘，矩阵相乘得到条带噪声，元素相乘得到斑点噪声
x_stripe=255*xin;
x_stripe=uint8(x_stripe);



% map = gray(256); % 256*3 数组中的每一行包含一种特定颜色的红、绿、蓝强度。强度介于 [0,1] 范围内，

% subplot(121)
% % imshow(255*xorig,map)
% imshow(uint8(255*xorig))
% title('original image')
% subplot(122)
% % imshow(x_stripe,map)
% imshow(x_stripe)
% title('noise contaminated image')