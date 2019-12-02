% covers the ref import
% set to zero or load existing values:
if ~exist(sprintf('%s\\%s\\ref_mean.mat',pwd,Dfolder),'file')
    ref_mean=zeros(nx,ny,numberofsets);   
else
    load(sprintf('%s\\%s\\ref_mean',pwd,Dfolder),'ref_mean');
end

disp(sprintf('skptrefimport: reference block nr _ out of %d processing',reftime/timestep));
for j = 1:reftime/timestep  %iterates over the blocks of data read at once
    disp([j])
    start = ((j-1)*timestep*freq)+1;
    stop = j*timestep*freq;
    ref = fread(fidref,timestep*freq*nx*ny,'uint16');
    ref = reshape(ref,nx,ny,timestep*freq);
    ref_mean(:,:,i) = ref_mean(:,:,i) + sum(ref,3);
end
ref_mean(:,:,i)=ref_mean(:,:,i)/(freq*reftime);
save(sprintf('%s\\%s\\ref_mean',pwd,Dfolder),'ref_mean');
clear ref
clear ref_mean
