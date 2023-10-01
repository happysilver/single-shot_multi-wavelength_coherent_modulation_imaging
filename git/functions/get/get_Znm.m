function Znm=get_Znm(zernike_num,N,m)

Znm = zeros(N,N,zernike_num);
for i= 1:zernike_num
%     Zt=get_Zernike(i+1,m);  %不要第一项
    mt=1.5*m;
    Zt=get_Zernike(i+1,mt).*get_Sig(m/2,mt,0,1);  %去掉外边缘看看
    Znm(:,:,i)=padarray(Zt,[(N-mt)/2,(N-mt)/2],0,'both');
%     figure,imshow(Znm(:,:,i),[])
end

end

function Znm= get_Zernike(j,res)
% 参数j：  Zernike多项式的序号
% 参数res：Zernike多项式的分辨率
    x           = linspace(-1,1,res);
    [x,y]       = meshgrid(x,x);
    [theta,rho] = cart2pol(x,y);% 由(x,y)换算(r,theta)
 
%     [n,m]       = Noll_to_Zernike(j); % 调用函数，返回n和m的值
    [n,m]     = OSA_to_Zernike(j-1);
%     [n,m]     = OSA_to_Zernike(j);  %不要第一项
    if m == 0
       deltam = 1;
    else
       deltam = 0;
    end
    Norm    = sqrt(2*(n+1)/(1+deltam));% 归一化因子 
    Rnm_rho = zeros(size(rho));        % 初始化径向多项式
    for s = 0:(n-abs(m))/2
        Rnm_rho = Rnm_rho+(-1)^s.*prod(1:(n-s))*rho.^(n-2*s)/(prod(1:s)*...
        prod(1:((n+abs(m))/2-s))*prod(1:((n-abs(m))/2-s))); % 径向多项式
    end
    if m < 0
       Znm = -Norm.*Rnm_rho.*sin(m.*theta); % m<0时候的zernike多项式
    else
       Znm = Norm.*Rnm_rho.*cos(m.*theta);  % m>0时候的zernike多项式
    end
    Znm = Znm.*(rho<=1); % mask = (rho<=1)，只保留单位圆内的数据
 end

function [n,m] = Noll_to_Zernike(j)
% 序号j从1开始；j = 1 对于的是piston模式
    if(j<1) 
        error('Noll indices start at 1.');
    end
    n  = 0;
    m  = 0;
    j1 = j-1;
    while (j1 > n)
        n  = n + 1;
        j1 = j1 - n;
        m  = (-1)^j * (mod(n,2) + 2*floor((j1+mod(n+1,2))/2));
    end
end

function [n,m] = OSA_to_Zernike(j)
 % 注意 这里的j是从0开始，j = 0对应的是piston模式，
 % 如果要从开始则传入的参数为j-1，此时可以与Noll_to_Zernike匹配
     n = floor(sqrt(2*j+1)+0.5)-1;
     m = 2*j-n*(n+2);
end