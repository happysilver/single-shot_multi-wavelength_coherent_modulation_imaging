close all;clear all;clc;

mkdir('..\..\dataset_new\TrainA');
mkdir('..\..\dataset_new\TrainB');
mkdir('..\..\dataset_new\TrainC');
addpath(genpath('..\..\dataset_new'));
addpath(genpath('..\..\functions'));

load('..\..\data_23-Mar-2023.mat','img_mono','object','img_poly','mask','delta','d','spect');

[size_m,size_n,num]=size(object);

all_num=1000;
[object_n,sup]=get_object(size_m,size_m/4,all_num-num,rand(8,8,all_num-num));
[img_poly_n,img_mono_n]=get_diffraction_img(object_n.*mask,delta,d,spect);
save('..\..\data_23-Mar-2023_extra.mat','img_mono_n','object_n','img_poly_n');

for nn=1:all_num
    filename1=['..\..\dataset_new\TrainA\A00'  num2str(nn) '.png'];
    filename2=['..\..\dataset_new\TrainB\B00'  num2str(nn) '.png'];
    filename3=['..\..\dataset_new\TrainC\C00'  num2str(nn) '.png'];
    if(nn<=num)
        imwrite(rescale(img_mono(:,:,nn)),filename1);
        imwrite(rescale(object(:,:,nn)),filename2);
        imwrite(rescale(img_poly(:,:,nn)),filename3);
    else
        imwrite(rescale(img_mono_n(:,:,nn-100)),filename1);
        imwrite(rescale(object_n(:,:,nn-100)),filename2);
        imwrite(rescale(img_poly_n(:,:,nn-100)),filename3);
    end
end