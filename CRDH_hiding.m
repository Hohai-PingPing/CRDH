clc
clear
addpath('./IWT/');
%========================自动生成秘密信息序列=========================
length=100000000; %定义秘密信息的二进制序列长度
payload=randi([0,1],length,1);
secret=payload';
%====================================================================
Carrier=imread('./imgs/1.Lena.pgm'); 
imwrite(uint8(Carrier),'./Carrier.png');
tic
[UN,UM,~]=size(Carrier);
Carrier_up=imresize(Carrier, [2*UN,2*UM],'nearest');
upsampling_time=toc;
disp(['上采样时间：',num2str(upsampling_time)]);

imwrite(uint8(Carrier_up),'./Carrier_up.png');
%==============================定义参数==============================
%====================================================================
[m,n,~]=size(Carrier_up);
mmmmm=m*n*(3/4);
%BPP=2;
%secret_length=floor(BPP*(m*n));
secret_length=100000;
secret=secret(1,1:secret_length);
save('secret.mat',"secret");
Carrier_up=double(Carrier_up);

tic
[A,H,V,D] = IWT_S(Carrier_up);
iwt_time=toc;
disp(['整数小波变换时间：',num2str(iwt_time)]);

tic
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
codebook_time=toc;
disp(['构建码本时间：',num2str(codebook_time)]);

tic
A_reshape=reshape(A,[1,(m/2)*(n/2)]);
for x=1:(m/2)*(n/2)
    if A_reshape(1,x)<=127
        symbolic(1,x)=-1;
    else
        symbolic(1,x)=1;
    end
end
symbolic=[symbolic,symbolic,symbolic];

bitBand=[reshape(r,[1,(m/2)*(n/2)]),reshape(b,[1,(m/2)*(n/2)]),reshape(a,[1,(m/2)*(n/2)])];

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

%位带bitband的前kk0个系数用于嵌入辅助信息，
%其中前kk个系数用于嵌入秘密信息的总长度，
%第kk+1到第kk0个系数用于嵌入a,b,r。
% kk0=24;
maxCapacity=sum(bitBand(1,kkbit+1:m*n*(3/4))); %该幅图像能够嵌入的最大容量
idealBPP=maxCapacity/(m*n);
disp(['该幅图像能够达到的最大容量为',num2str(maxCapacity)]);
weiweiwei=max(size(dec2bin(maxCapacity)));
capacityInformation=dec2bin(secret_length,idealbits);%代表容量总长度的秘密信息

if secret_length<=maxCapacity
    disp(['该幅图像嵌入了长度为',num2str(secret_length),'的秘密信息']);
else
    disp(['该幅图像的最大容量为',num2str(maxCapacity),'，不能将长度为',num2str(secret_length),'的秘密信息全部嵌入']);
end

label=0;
b=1;
for y=1:7
    for g=kkbit+1:m*n*(3/4)
        if y<=bitBand(1,g)
            embed(g,y)=secret(1,b);
            b=b+1;
        else
            continue;
        end
        if b-1>=secret_length
            label=1;
            break;
        end
    end

    if label==1
        break;
    end
end

actualEmbeddedCapacity=b-1;
disp(['实际嵌入容量为',num2str(actualEmbeddedCapacity)]);

jilun=y;
aitemp=g;

if y>1
    templength=m*n*(3/4);
else
    templength=g;
end

for f=kkbit+1:templength
    if f<=g
        temp000=min(jilun,bitBand(1,f));
        if temp000==0
            Eembed(f,1)=0;
        else
        Eembed(f,1)=symbolic(1,f)*bin2dec(num2str(fliplr(embed(f,1:temp000))));
        end
    else
        temp000=min(jilun-1,bitBand(1,f));
        if temp000==0
            Eembed(f,1)=0;
        else
        Eembed(f,1)=symbolic(1,f)*bin2dec(num2str(fliplr(embed(f,1:temp000))));
        end
    end
end

AI=capacityInformation;
label=0;
b=1;
for y=1:sbits
    for g=1:kkbit
        if y<=bitBand(1,g)
            aaa=AI(1,b);
            embed(g,y)=str2num(AI(1,b));
            b=b+1;
        else
            continue;
        end

        if b-1>=idealbits
            label=1;
            break;
        end
    end
    if label==1
        break;
    end
end
jilun0=y;
aitemp0=g;
if y>1
    templength=kkbit;
else
    templength=g;
end
for f=1:kkbit
    if f<=aitemp0
        temp000=min(jilun0,bitBand(1,f));
        if temp000==0
            Eembed(f,1)=0;
        else
        Eembed(f,1)=symbolic(1,f)*bin2dec(num2str(fliplr(embed(f,1:temp000))));
        end
    else
        temp000=min(jilun0-1,bitBand(1,f));
        if temp000==0
            Eembed(f,1)=0;
        else
            Eembed(f,1)=symbolic(1,f)*bin2dec(num2str(fliplr(embed(f,1:temp000))));
        end
    end
end

coefficient_length=max(size(Eembed));

save('Eembed.mat',"Eembed");

Dd=zeros([1,m*n*(1/4)]);
Vv=zeros([1,m*n*(1/4)]);
Hh=zeros([1,m*n*(1/4)]);

zonglength=coefficient_length;
for oo=1:3
    Dd(1,1:min(zonglength,(m*n*(1/4))))=Eembed(1:min(zonglength,(m*n*(1/4))),1);
    if zonglength>(m*n*(1/4))
        Vv(1,1:min(zonglength,m*n*(2/4))-(m*n*(1/4)+1)+1)=Eembed(m*n*(1/4)+1:min(zonglength,m*n*(2/4)),1);
        if zonglength>(m*n*(2/4))
            Hh(1,1:min(zonglength,m*n*(3/4))-(m*n*(2/4)+1)+1)=Eembed(m*n*(2/4)+1:min(zonglength,m*n*(3/4)),1);
        else
            break;
        end
    else
        break;
    end
end

DD=reshape(Dd,[m/2,n/2]);
VV=reshape(Vv,[m/2,n/2]);
HH=reshape(Hh,[m/2,n/2]);

stego=RIWT_S(A,HH,VV,DD);

embedding_time=toc;
disp(['信息嵌入时间：',num2str(embedding_time)]);

total_time=upsampling_time+iwt_time+codebook_time+embedding_time;
disp(['秘密信息隐藏的总时间：',num2str(total_time)]);

a1=upsampling_time/total_time;
a2=iwt_time/total_time;
a3=codebook_time/total_time;
a4=embedding_time/total_time;

PSNR=psnr(uint8(Carrier_up),uint8(stego));
disp(['Stego-image的PSNR=',num2str(PSNR),'dB']);
save('stego.mat',"stego");

for i=1:m
    for j=1:n
        if (stego(i,j)>255)||(stego(i,j)<0)
            disp(['第',num2str(i),',',num2str(j),'个像素值溢出']);
        end
    end
end
imwrite(uint8(stego),'./stegotest.png');

figure
subplot(121);imshow(uint8(Carrier_up)),title('上采样图像');
subplot(122);imshow(uint8(stego)),title('伪装图像');
















