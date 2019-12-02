%% Initialisation of variables

if ~exist(STDfolder,'dir')
    mkdir(STDfolder);
end

if ~exist(sprintf('%s\\%s\\D_std.mat',pwd,STDfolder),'file')
    D_std=zeros(nx,ny,numberofsets);      % will contain the STD value for each meshpoint
else
    load(sprintf('%s\\%s\\D_std',pwd,STDfolder),'D_std')
end
load(sprintf('%s\\%s\\D_mean',pwd,Dfolder),'D_mean');
D=zeros(nx,ny,freq*timestep);           % mixing scalar
disp(sprintf('STD set %d frames _ to _ out of %d',i,singletime*freq));
%% STD
%status = fseek(fidsingle,0,'bof'); %set the pointer back to the beginning of file
for j=1:numberofblocks  %iterates the blocks of data read in at once
    start = ((j-1)*timestep*freq)+1;    %only for display reasons
    stop = j*timestep*freq;             %
    disp([start stop]);                 %
    load(sprintf('%s\\%s\\D_set%02d_block%02d',pwd,Dfolder,i,j),'D')
    
    for l=1:timestep*freq       % iterates the timeframes of one block of data
        D_std(:,:,i)=(D(:,:,l)-D_mean(:,:,i)).^2+D_std(:,:,i);
    end
end

D_std(:,:,i)=(sqrt(D_std(:,:,i)/(singletime*freq-1))).*mask;
save(sprintf('%s\\%s\\D_std',pwd,STDfolder),'D_std');

% %% plot
% createfigurestd(D_std(:,:,i))
% saveas(gcf,sprintf('%s\\%s\\STD_%d',pwd,STDfolder,i),'jpg');
%% clean up
close gcf
clear D
clear D_std
clear start stop
clear D_mean


