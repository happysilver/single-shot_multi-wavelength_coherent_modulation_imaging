%% This function is to get circle region of the matrix.
function Sig = get_Sig( R, N, out, in )
    x=-55:110/(N-1):55;
    y=x;
    Sig=out*ones(N,N);
    if R<0
        Sig=out.*ones(N,N);
        Sig(-R+1:N+R,-R+1:N+R)=in;
    else
        R=55*2*R/N;
        for i=1:N
            for j=1:N
                xx=x(i);
                yy=y(j);
                p=xx^2+yy^2;
                if p<=R^2
                    Sig(i,j)=in;
                end
            end
        end
    end
end