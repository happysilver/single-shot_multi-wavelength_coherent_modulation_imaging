function [object,sup]=get_object(N,m,object_num,sample)

object=zeros(N,N,object_num);
sup=get_Sig(m/2,N,0,1);
for i=1:1:object_num
    object0=imresize(sample(:,:,i),[m m]);
    object1=padarray(object0,[m/2*3 m/2*3],0,'both');
    object(:,:,i)=object1.*sup;
end

end