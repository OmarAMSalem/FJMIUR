function [selected_fs]=FJMIUR(data)
[newdata] = ISUR( data );
[selected_fs]=FJMI(newdata);
end

function [best_fs]=FJMI(data)
data(:,end) = grp2idx(data(:,end));
NC = max(data(:,end));
NF = size(data,2) - 1;
NS=size(data,1);
I_Cx = zeros(1,NF);
H_x = sparse(1,NF);
I_Cx = sparse(1,NF);
H_xx = sparse(NF,NF);
I_xx = sparse(NF,NF);
SN=sparse(NS,NS);
H_SN=0;
newdata=[];

[relc]= decision(data(:,end));
H_C =  -(sum(log2(sum(relc)/length(relc))))/length(relc);

max_MI=0;firstFeature=1;
for i=1:NF
[relx]= relation(data(:,i));
 I_Cx(i)= F_MI(relx,relc);

    if I_Cx(i)> max_MI
        max_MI=I_Cx(i);
        firstFeature=i;
    end
end

best_fs(1)=firstFeature;
best_val(1)=max_MI;
SN=relation(data(:,firstFeature));

selected=zeros(1,NF);
selected(best_fs(1))=1;

for n=2:size(data,2)-1
    max_IcSNx=0;
    bestFeature=0;
    for i=1:NF
       if selected(i) continue;end;
       [relxn]= relation(data(:,i));  
    I_SNxc=iF_MI(relxn,SN,relc);
     I_Cx=F_MI(relxn,relc);
     I_CSN=F_MI(SN,relc);
     JRes=I_Cx+I_CSN-I_SNxc;
       if JRes>max_IcSNx
           max_IcSNx=JRes;
           bestFeature=i;
       end
       
    end

    if max_IcSNx <= best_val(n-1) 
         break; 
    end

    best_fs(n)=bestFeature;
    SN=min(SN,relation(data(:,bestFeature)));
    selected(bestFeature)=1;
    best_val(n)=max_IcSNx;
end
end
 


function [rel]= relation(x)
for i=1:length(x)
    for j=1:length(x)
        rel(j,i)=exp(-abs(x(i)-x(j)));
    end
end
end


function [rel]= decision(x)
for i=1:length(x)
    for j=1:length(x)
        if x(i)==x(j)
         rel(j,i)= 1;
        else rel(j,i)= 0;
        end
  
    end
end
end

function [fMI]= F_MI(relx,rely)
    H_x =  -(sum(log2(sum(relx)/length(relx))))/length(relx);
    H_y = -(sum(log2(sum(rely)/length(rely))))/length(rely);
    relxy = min(rely,relx); 
    H_xy=  -(sum(log2(sum(relxy)/length(relxy))))/length(relxy);
    fMI = (H_x + H_y - H_xy);
end

function [IfMI]= iF_MI(relx,rely,relz)
    H_x =  -(sum(log2(sum(relx)/length(relx))))/length(relx);
    H_y = -(sum(log2(sum(rely)/length(rely))))/length(rely);
    H_z = -(sum(log2(sum(relz)/length(relz))))/length(relz);
    
    relxy = min(rely,relx); 
    H_xy=  -(sum(log2(sum(relxy)/length(relxy))))/length(relxy);
    relxz = min(relz,relx); 
    H_xz=  -(sum(log2(sum(relxz)/length(relxz))))/length(relxz);
    relyz = min(relz,rely); 
    H_yz=  -(sum(log2(sum(relyz)/length(relyz))))/length(relyz);
    
    relxyz=min(relxy,relz);
    H_xyz=  -(sum(log2(sum(relxyz)/length(relxyz))))/length(relxyz);
    
    IfMI = (H_x + H_y + H_z + H_xyz)-(H_xy+H_xz+H_yz);
end

function [ newdata ] = ISUR( data )
class=data(:,end);
uniclass=unique(class); 
impdata=[data(:,1:end-1) class];
classratio = histc(class,uniclass);
[aval,bindx]=min(classratio); 	
mindata=impdata(find(class == bindx),:);
maxdata=impdata(find(class ~= bindx),:);
[sel_indx] = distfunc(mindata(:,1:end-1),maxdata(:,1:end-1));
mndata=data(find(class == bindx),:);
mxdata=data(find(class ~= bindx),:);
newdata=[mndata;mxdata(sel_indx,:)];
end

function [selected] = distfunc(X,Y)
Xsize=size(X,1);
D = pdist2(X,Y);
D2=D';
selected=zeros(1,Xsize);

for i=1:Xsize
    [a,b]=min(D2(:,i));
    selected(i)=b;
D2(b,:)=inf;

end
end
