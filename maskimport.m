mask=importdata('mask.asc');
%subplot(1,2,1)
%surf(mask)

for i=1:128
    for j=1:64
        if mask(i,j)==0
            mask(i,j)=NaN;
        end
    end
end
clear i j;
%subplot(1,2,2)
%surf(mask)