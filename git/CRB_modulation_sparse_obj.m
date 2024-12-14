function CRB_modulation_sparse_obj()
close all;clear;clc

addpath(genpath('functions'));

d=225e-3;
delta=6.5e-6;
spect=[405 488 549 610 638]*1e-9;

N = 256; % size of the diffraction pattern image
m = 64; % size of the sample image

x_num=3;
theta0=[2*pi/40,pi/7,0.5];
[xx,yy]=meshgrid(-N/2:(N+1)/N:N/2);
A=exp(-(xx.^2+yy.^2)/(m/2)^2);

x0=rand(16,16);

tic
% options=optimoptions(@fmincon,);
x=fmincon(@least_func,x0,[],[],[],[],zeros(size(x0)),ones(size(x0)),[]);
toc

mask=exp(2i*pi*imresize(x,[N,N]));
obj=make_obj(theta0);
I1=get_poly_img(obj);
I2=get_poly_img(obj.*mask);

figure,imshow(angle(obj),[])
figure,imshow(abs(obj),[])

figure,
subplot(1,2,1),imshow(log(1+I1),[])
subplot(1,2,2),imshow(log(1+I2),[])
colormap("parula")


save_file=['mask.mat'];
save(save_file,"mask","theta0","obj","I1","I2","delta","d","spect");


disp("end")

    function CRB=least_func(x)
        % x -> mask
        mask=exp(2i*pi*imresize(x,[N,N]));

        delta_theta=1e-3;
        reg_para=1e-4;

        H=zeros(N*N,x_num);
        I_c=get_poly_img(make_obj(theta0).*mask);
        for repeat_least_func=1:x_num
            theta_p=theta0;
            theta_p(repeat_least_func)=theta_p(repeat_least_func)+delta_theta;
            theta_n=theta0;
            theta_n(repeat_least_func)=theta_n(repeat_least_func)-delta_theta;
            I_p=get_poly_img(make_obj(theta_p).*mask);
            I_n=get_poly_img(make_obj(theta_n).*mask);
            H_t=(I_c+reg_para).^(-1/2).*abs(I_p-I_n)/(2*delta_theta);
            H(:,repeat_least_func)=H_t(:);
        end
        Fisher_matrix=H'*H;
        CRB=sum(Fisher_matrix(:));
    end

    function I_out=get_poly_img(obj)
        I_out=zeros(N,N);
        for repeat_get_poly_img=1:1:length(spect)
            I_out=I_out+abs(fresnel_advance(obj, delta, delta, d, spect(repeat_get_poly_img))).^2;
        end
    end

    function y=make_obj(theta)
        omega=theta(1);
        bias=theta(2);
        amplitude=theta(3);
        %         rotate_angle=theta(4);

        x_t=1:1:m;
        grating=amplitude*sin(omega.*x_t+bias);
        obj_t=repmat(grating,[m,1]);
        obj_n=padarray(obj_t,[(N-m)/2,(N-m)/2],0);
        %         obj_n=imrotate(obj_n,rotate_angle/pi*180,"bilinear","crop");
        y=A.*exp(1i*obj_n);
    end
end