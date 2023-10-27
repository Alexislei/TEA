

function off=repair_c(off,ch)

    %处理empty space
    [m,n]=size(ch);
    zeronum_ch=sum(sum(ch==0));
    zeronum_off=sum(sum(off==0));
    
    %处理empty位置
    missItem=setdiff(ch,off); % off 中未出现的元素
    multiItem=find_dup(off); % off中重复的元素，包含empty
    
    if zeronum_off<zeronum_ch
        missItem(missItem==0)=[];
        for i=1:(zeronum_ch-zeronum_off)
            missItem=[missItem;0]; %补全消失的empty位置
        end
        randp=randperm(size(missItem,1));
        missItem=missItem(randp');
    end    
    
%     if zeronum_off<zeronum_ch
%         size_mI=size(missItem,1);
%         zeronum_mI=sum(sum(missItem==0));
%         for i=1:(zeronum_ch-zeronum_off-zeronum_mI)
%             missItem(size_mI+i,:)=0;
%         end
% %         missItem(end+1:end+(zeronum_ch-zeronum_off),:)=0; %补全消失的empty位置
%         randp=randperm(size(missItem,1))';
%         missItem=missItem(randp');
%     end

    if ~isempty(missItem)
        mask=zeros(m,n);
        for i=1:size(multiItem,1)
            multi_ma=reshape((off==multiItem(i,:)),[1,numel(off)])';
            idx=find(multi_ma==1);
            if multiItem(i,:)==0
                if size(idx,1)>zeronum_ch || size(idx,1)==zeronum_ch
                    multi_ma(idx(randperm(size(idx,1),zeronum_ch),:),:)=0; %重复元素只保留zeronum_ch个，保留的那个元素置为0
                else
                    multi_ma=zeros(m,n);
                end
            else
                multi_ma(idx(randperm(size(idx,1),1),:),:)=0; %重复元素只保留一个，保留的那个元素置为0
            end 
            mask=mask+reshape(multi_ma,[m,n]);
        end
        mask_val=mask;
        mask_val(mask_val==1)=missItem;
        off=(1-mask).*off+mask.*mask_val;
    end
end

function dup_item=find_dup(a)
    b = unique(a);
%     b(b==0)=[]; 
%     b(b==-1)=[];
    nn = length(b);
    pos_dup = [];
    dup_item=[];
    for i=1:nn
        pos=[];
        num_b = b(i,1);
        [line1, row1] = find(a==num_b);
        if size(line1,1)>1
            pos(:,1)= line1;
            pos(:,2)= row1;
            dup_item=[dup_item;num_b];
        end
%         pos_dup=[pos_dup;pos];
    end
end
