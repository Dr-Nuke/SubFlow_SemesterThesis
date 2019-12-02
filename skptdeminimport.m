% covers the demin import
% set to zero or load existing values:
if ~exist(sprintf('%s\\%s\\demin_mean.mat',pwd,Dfolder),'file')
    demin_mean=zeros(nx,ny,numberofsets);   
else
    load(sprintf('%s\\%s\\demin_mean',pwd,Dfolder),'demin_mean');
end

disp(sprintf('skptdeminimport: demin block nr _ out of %i processing',demintime/timestep));
for j = 1:demintime/timestep  %iterates over the blocks of data read at once
    disp([j])
    start = ((j-1)*timestep*freq)+1;
    stop = j*timestep*freq;
    demin = fread(fiddemin,timestep*freq*nx*ny,'uint16');
    demin = reshape(demin,nx,ny,timestep*freq);
    demin_mean(:,:,i) = demin_mean(:,:,i)+sum(demin,3);
end
demin_mean(:,:,i)=demin_mean(:,:,i)/(freq*demintime);
save(sprintf('%s\\%s\\demin_mean',pwd,Dfolder),'demin_mean');
clear demin
clear demin_mean


