function zernike_obj_modulated_CDI

close all;clear;clc

addpath(genpath('functions'));

d=225e-3;
delta=6.5e-6;
spect=[405,488,638]*1e-9;
lambda_c=mean(spect);
spect_num=length(spect);
spect_I=normpdf(spect./lambda_c,1,0.3);
spect_I=spect_I./max(spect_I(:));

diffraction_option='mono';

N = 1024; % size of the diffraction pattern image
m = 512; % size of the sample image

zernike_num = 21-1;
Znm=get_Znm(zernike_num,N,m);
sup=SUP;
obj_theta=rand(1,zernike_num)*0.5+0.25;

% add block
block_area=ones(N,N);
% [xt,yt]=meshgrid(1:1:N);
% block_area_r=N/10;
% block_area=((xt-N/2).^2+(yt-N/2).^2)>block_area_r.^2;


% mask
flag_mask=1;
load('mask.mat',"mask");
mask=imresize(mask,[m,m]);
[mt,nt]=size(mask);
mask=padarray(mask,[(N-mt)/2,(N-nt)/2],0);
figure,imshow(mask,[]),colormap("parula")

obj=make_obj(obj_theta,Znm);
x_ori(1,1:zernike_num)=obj_theta;
x_ori(2,1:spect_num)=spect;
x_ori(3,1:spect_num)=spect_I;
I0=get_diffraction_img(x_ori,'mono');
I1=get_diffraction_img(x_ori,diffraction_option);


adjust dynamic range
DR=65535;
I1=uint16(I1./max(I1(:))*DR);

I1=double(I1);
I1=I1./max(I1(:));

x0(1,1:zernike_num)=rand(1,zernike_num);
x0(2,1:spect_num)=spect_0;
x0(3,1:spect_num)=spect_I0;
least_func_alpha=1;

tic
x=x0;
x=ga(@least_func,x0,[],[],[],[],zeros(size(x0)),ones(size(x0)),[]);
toc

obj_out=make_obj(x(1,:),Znm).*sup;

I_out=get_diffraction_img(x,diffraction_option);
if(flag_mask)
    obj_show=obj./mask;
else
    obj_show=obj;
end

spect_out=x(2,1:spect_num);
spect_I_out=x(3,1:spect_num);

obj_range=N/2-m/2:1:N/2+m/2;

figure,
subplot(3,3,1),imshow(angle(obj_show(obj_range,obj_range)),[]),title('object')
subplot(3,3,2),imshow(I1,[]),title('poly image')
% subplot(3,3,3),imshow(I1-I_noise,[]),title('no noise image')
subplot(3,3,3),imshow(I0,[]),title('mono image')
subplot(3,3,4),imshow(angle(obj_out(obj_range,obj_range)),[]),title('restored object')
subplot(3,3,5),imshow(I_out,[]),title('restored image')
subplot(3,3,7),plot(1:zernike_num,x_ori(1,1:zernike_num),1:zernike_num,x(1,1:zernike_num)),title(['SNR = ',num2str(round(SNR))]),legend('original x','restored x'),ylim([0,1])
subplot(3,3,8),plot(abs(x(1,1:zernike_num)-x_ori(1,1:zernike_num))./x_ori(1,1:zernike_num)*100,'*'),ylabel("relative error(%)")
subplot(3,3,9),hold on,plot(spect/1e-6,spect_I,'o'),plot(spect_0/1e-6,spect_I0,'*'),plot(spect_out/1e-6,spect_I_out,'*'),legend('original spect','initial guess','restored spect'),hold off
colormap("parula")

disp('end')

    function I_out=get_diffraction_img(x,diffraction_option)
        if(flag_mask)
            obj=make_obj(x(1,:),Znm).*mask;
        else
            obj=make_obj(x(1,:),Znm);
        end
        if(diffraction_option=='mono')
            I_out=abs(fresnel_advance(obj.*sup, delta, delta, d, lambda_c)).^2;
        elseif(diffraction_option=='poly')
            I_out=zeros(N,N);
            for get_diffraction_img_index=1:spect_num
                I_out=I_out+x(3,get_diffraction_img_index).*abs(fresnel_advance(obj.*sup, delta, delta, d, x(2,get_diffraction_img_index))).^2;
            end
        end
        I_out=I_out.*block_area;
    end

    function y=least_func(x)
        I_temp=get_diffraction_img(x,diffraction_option);
        Regular_term=@(x)(norm(gradient(x),2));
        log_likelihood=I_temp-I1.*log(handle_zero(I_temp));
        y=sum(log_likelihood(:))+least_func_alpha*Regular_term(I_temp);
%         y=norm(I_temp-I1);
    end

    function [y,grad]=fun(x)
        y=get_diffraction_img(x,diffraction_option);

        grad=zeros(1,zernike_num);
        delta_x=1e-3;
        for repeat_fun=1:zernike_num
            x_p=x;
            x_p(repeat_fun)=x_p(repeat_fun)+delta_x;
            x_n=x;
            x_n(repeat_fun)=x_n(repeat_fun)-delta_x;
            grad(repeat_fun)=(get_diffraction_img(x_p,diffraction_option)-get_diffraction_img(x_n,diffraction_option))/(2*delta_x);
        end        
    end

end

function y=make_obj(theta,Znm)

[N,~,zernike_num]=size(Znm);
y=zeros(N,N);
for i=1:zernike_num
    y=y+Znm(:,:,i).*theta(i)/zernike_num;
end
y=exp(1i*pi.*y);

end

function I_out=handle_zero(I_in)
    I_out=I_in;
    I_out(I_out==0)=1;
end