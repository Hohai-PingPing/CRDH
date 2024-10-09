clc
clear
originalSecret=load('./secret.mat');
originalSecret=originalSecret.secret;
originalLength=max(size(originalSecret));

Carrier=imread('./Carrier.png');
stego=imread('./stegotest.png');

tic
[m,n,~]=size(stego);
mmmm=m*n*(3/4);
[A,H,V,D] = IWT_S(stego);

PSNR_cover=psnr(uint8(Carrier),uint8(A));

for i=1:m/2
    for j=1:n/2
        if A(i,j)<=127
            a(i,j)=floor((log2(257-A(i,j)))-1);
            b(i,j)=a(i,j);
            r(i,j)=floor(log2(A(i,j)+1));
        else
            a(i,j)=floor((log2(A(i,j)+2))-1);
            b(i,j)=a(i,j);
            r(i,j)=floor(log2(256-A(i,j)));
        end
    end
end
bitBand=[reshape(r,[1,(m/2)*(n/2)]),reshape(b,[1,(m/2)*(n/2)]),reshape(a,[1,(m/2)*(n/2)])];

informationBand=[reshape(D,[1,m*n*(1/4)]),reshape(V,[1,m*n*(1/4)]),reshape(H,[1,m*n*(1/4)])];
informationBand=abs(informationBand);

idealCapacity=sum(bitBand(1,1:m*n*(3/4)));
idealbits=max(size(dec2bin(idealCapacity)));

sbits=1;
temp1=0;
kkbit=0;
for k=1:(m*n*(3/4))
    if temp1>=idealbits
        kk=k-1;
        break;
    else
        temp1=temp1+min(sbits,bitBand(1,k));
    end
end
kkbit=kk;

%===============提取辅助信息================
yiulunLength=kkbit;
AI='';
label=0;
for kkk=1:sbits
    for k=1:yiulunLength
        length=max(size(AI));
        if length>=idealbits
            label=1;
            break;
        end
            aaa=dec2bin(informationBand(1,k),min(sbits,bitBand(1,k)));
            if kkk<=bitBand(1,k)
                AI=strcat(AI,aaa(1,min(sbits,bitBand(1,k))+1-kkk));
            else
                continue;
            end

    end %for k=22:m*n*(3/4)
    if label==1
        break;
    end
end

secret_length=bin2dec(AI);
%=============提取秘密信息=================
yiulunLength=m*n*(3/4);
secret='';
label=0;
jjj=1;
for kkk=1:7
    for k=kkbit+1:yiulunLength
        length=max(size(secret));
        if length>=secret_length
            label=1;
            break;
        end

            aaa=dec2bin(informationBand(1,k),min(8,bitBand(1,k)));
            if kkk<=bitBand(1,k)
                secret=strcat(secret,aaa(1,min(8,bitBand(1,k))+1-kkk));
                jjj=jjj+1;
            else
                continue;
            end
    end %for k=22:m*n*(3/4)
    if label==1
        break;
    end
end

extraction_time=toc;
disp(['提取时间：',num2str(extraction_time)]);

total=0;
if originalLength==length
    for h=1:length
        secretMatrix(1,h)=str2num(secret(1,h));
        if originalSecret(1,h)==secretMatrix(1,h)
            total=total+1;
        else
            disp(['第',num2str(h),'位出错']);
        end   
    end
    accuracy=total/length;
    if accuracy==1
        disp('秘密信息100%恢复');
    else
        disp(['秘密信息恢复出错，其精确度为',num2str(accuracy)]);
    end
else
    disp(['秘密信息恢复严重出错，其长度为',num2str(length),',而理论长度为',num2str(originalLength)]);
end
























