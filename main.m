clear;
clc;
close all;

hight=100;

%Importing original image
im=imread('Football-grass2.jpg');
im=rgb2gray(im);
im=im2double(im);
[m0,n0]=size(im);
scaling=hight/m0;

%Smoothing and resizing 
h=fspecial('average',1);
im2=imfilter(im,h,'replicate');
im2=imresize(im2,scaling);
[m,n]=size(im2);
N=m*n;

%Calculating distance between pixels
disp('Calculating distance between pixels.')
x=0:m-1;
y=0:n-1;
xpos=repmat(x',1,n);
ypos=repmat(y,m,1);
xpos=xpos(:);
ypos=ypos(:);
d=sqrt((repmat(xpos,1,N)-repmat(xpos',N,1)).^2+...
    (repmat(ypos,1,N)-repmat(ypos',N,1)).^2);
clear xpos y_pos x y

%Calculating difference in intensity between pixels
disp('Calculating difference in intensity.')
diff_im=im2(:)*ones(1,N)-ones(1,N)'*im2(:)';

%Find pixels closer to each other than r
r=16;
valid=sparse(d<r);

%Calculating W
disp('Calculating W.')
sigma1=1/10;
sigmad=1.3;
W=zeros(N);
W(valid)=exp(-diff_im(valid).^2/sigma1^2-d(valid)/sigmad^2);
clear d diff_im valid D
W=sparse(W);

%Calculating A in eigenvalue problem
D=diag(sum(W,2));
D2=sqrt(inv(D));
A=D2*(D-W)*D2;

%Solving eigenvalue problem
disp('Solving eigenvalue problem.')
opts.IsFunctionSymmetric=1;
[Z,E]=eigs(A,2,'smallestabs',opts);%Assuming A is positiv defined
y=D2*Z(:,2);
clear Z D

%Finding some cutting points
B=50;
s=linspace(min(y)+eps,max(y),B);

%Calculating the best Ncut
disp('Calcualting the best Ncut.')
Ncut=zeros(1,B);
for i=1:B
    a=y<s(i);
    cutAB=sum(sum(W(a,~a)));
    assocAV=sum(sum(W(a,:)));
    assocBV=sum(sum(W(~a,:)));
    Ncut(i)=cutAB/assocAV+cutAB/assocBV;
end
[~,min_i]=min(Ncut);

%Making the new image and showing it
a=y<s(min_i);
im3=reshape(a,[m,n]);
subplot(1,3,1),imshow(im);
title('Original')
subplot(1,3,2),imshow(im2);
title('Scaled')
subplot(1,3,3),imshow(im3);
title('Segmentation')
disp('Done')