function sparse_obj_modulated_CDI

close all;clear;clc

addpath(genpath('functions'));

d=225e-3;
delta=6.5e-6;
spect=[405,488,638]*1e-9;
lambda_c=mean(spect);
spect_num=length(spect);
spect_I=normpdf(spect./lambda_c,1,0.3);
spect_I=spect_I./max(spect_I(:));

diffraction_option='poly';

N = 1024; % size of the diffraction pattern image
m = 512; % size of the sample image

%% add block
block_area=ones(N,N);
% [xt,yt]=meshgrid(1:1:N);
% block_area_r=N/10;
% block_area=((xt-N/2).^2+(yt-N/2).^2)>block_area_r.^2;

%% mask
flag_mask=1;
load('mask.mat',"mask");
mask=imresize(mask,[m,m]);
[mt,nt]=size(mask);
mask=padarray(mask,[(N-mt)/2,(N-nt)/2],0);
figure,imshow(angle(mask),[]),colormap("parula")

%% obj
% load('source.mat','ss');
A=1;

x_ori=[2*pi/80,pi/5];
obj=make_obj(x_ori);
I1=get_diffraction_img(x_ori,diffraction_option);

x_num=2;

tic
x=ga(@least_func,x_num,[],[],[],[],zeros(1,x_num),[pi/20,pi]);
toc

obj_out=make_obj(x);
obj_range=N/2-m/2:1:N/2+m/2;
I_out=get_diffraction_img(x,diffraction_option);
obj_angle=angle(obj);l0=obj_angle(N/2,(N-m)/2+1:(N+m)/2);
obj_out_angle=angle(obj_out);l1=obj_out_angle(N/2,(N-m)/2+1:(N+m)/2);

figure,imshow(angle(obj),[]),colormap("parula")
figure,imshow(I1,[]),colormap("turbo")

figure,
subplot(3,2,1),imshow(angle(obj(obj_range,obj_range)),[]),title('object')
subplot(3,2,2),imshow(angle(obj_out(obj_range,obj_range)),[]),title('restored object')
subplot(3,2,3),imshow(I1,[]),title('poly image')
subplot(3,2,4),imshow(I_out,[]),title('restored image')
subplot(3,1,3),hold on,plot(l0),plot(l1),legend('object','restored object')
colormap("parula")


disp('end')

    function I_out=get_diffraction_img(x,diffraction_option)
        if(flag_mask)
            obj_t=make_obj(x).*mask;
        else
            obj_t=make_obj(x);
        end
        if(diffraction_option=='mono')
            I_out=abs(fresnel_advance(obj_t, delta, delta, d, lambda_c)).^2;
        elseif(diffraction_option=='poly')
            I_out=zeros(N,N);
            for get_diffraction_img_index=1:spect_num
                I_out=I_out+spect_I(get_diffraction_img_index).*abs(fresnel_advance(obj_t, delta, delta, d, spect(get_diffraction_img_index))).^2;
            end
        end
        I_out=I_out.*block_area;
        I_out=DR_img(I_out);
    end

    function y=least_func(x)
        I_temp=get_diffraction_img(x,diffraction_option);
%         Regular_term=@(x)(norm(gradient(x),2));
%         log_likelihood=I_temp-I1.*log(handle_zero(I_temp));
%         y=sum(log_likelihood(:))+least_func_alpha*Regular_term(I_temp);
        y=norm(I_temp-I1);
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
    function y=make_obj(theta)
        omega=theta(1);
        bias=theta(2);
        % amplitude=theta(3);
%         rotate_angle=theta(4);

        x_t=1:1:m;
        % grating=amplitude*sin(omega.*x_t+bias);
        grating=sin(omega.*x_t+bias);
        obj_t=repmat(grating,[m,1]);
        obj_n=padarray(obj_t,[(N-m)/2,(N-m)/2],0);
%         obj_n=imrotate(obj_n,rotate_angle/pi*180,"bilinear","crop");
        y=A.*exp(1i*obj_n);
    end

end

function I_out=handle_zero(I_in)
    I_out=I_in;
    I_out(I_out==0)=1;
end

function I_out=DR_img(I_in)
%     DR=2^20-1;
%     I_out=(I_in-min(I_in(:)))./(max(I_in(:))-min(I_in(:))); %归一化到0-1
%     I_out=double(round(I_out*DR));
    
    I_out=I_in;
end