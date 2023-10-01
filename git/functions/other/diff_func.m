function diff_out=diff_func(I1,I2,option)
    I1=MatMap(I1,0,1);
    I2=MatMap(I2,0,1);
    if(option=="corr")
        diff_out=corr2(I1,I2);  
    elseif(option=="var")
        diff_temp=(I1-I2).^2;
        diff_out=sum(diff_temp(:));
    end
end