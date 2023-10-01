function [img_poly,img_mono]=get_diffraction_img(object,delta,d,spect)

[m,n,object_num]=size(object);
img_poly=zeros(size(object));
img_mono=zeros(size(object));
for i=1:object_num
    u0=object(:,:,i);
    % poly
    img_sum=zeros(m,n);
    img_temp=zeros(m,n);
    for lambda=spect
        img_temp=fresnel_advance(u0, delta, delta, d, lambda);
        img_temp=abs(img_temp);
        img_sum=img_sum+img_temp;
    end
    img_poly(:,:,i)=img_sum;

    %mono
    lambda=mean(spect);
    img_temp=fresnel_advance(u0, delta, delta, d, lambda);
    img_temp=abs(img_temp);
    img_mono(:,:,i)=img_temp;
end

end