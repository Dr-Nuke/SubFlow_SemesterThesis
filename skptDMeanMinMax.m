% creates D_min, D_mean, D_max

% memory is allocated if not already done.
% eventually previous generated data is loaded


%% initializing variables, load or create D_minmeanmax

if ~exist(sprintf('%s\\%s\\D_min.mat',pwd,Dfolder),'file')
    D_min=ones(nx,ny,numberofsets);   
    %map of min values for any meshpoint (sensor 1)
else
    load(sprintf('%s\\%s\\D_min',pwd,Dfolder),'D_min');
end
if ~exist(sprintf('%s\\%s\\D_mean.mat',pwd,Dfolder),'file')
    D_mean=zeros(nx,ny,numberofsets);   
    %map of mean values for any meshpoint (sensor 1)
else
    load(sprintf('%s\\%s\\D_mean',pwd,Dfolder),'D_mean');
end
if ~exist(sprintf('%s\\%s\\D_max.mat',pwd,Dfolder),'file')
    D_max=zeros(nx,ny,numberofsets);   
    %map of max values for any meshpoint (sensor 1)
else
    load(sprintf('%s\\%s\\D_max',pwd,Dfolder),'D_max');
end
load(sprintf('%s\\%s\\ref_mean',pwd,Dfolder),'ref_mean');
load(sprintf('%s\\%s\\demin_mean',pwd,Dfolder),'demin_mean');
D_meantemp = zeros(nx,ny);
D=zeros(nx,ny,freq*timestep);
onet=ones(nx,ny);

%% computations
disp(sprintf(...
    'skpDMeanMax: processing single set %d frame _ to _ out of %d',...
    i,freq*singletime));
status = fseek(fidsingle,0,'bof');  % set the pointer back to the
                                    % beginning of file
for j = 1:numberofblocks        % iterates the data blocks
    start = ((j-1)*timestep*freq)+1;
    stop = j*timestep*freq;
    disp([start stop]);
    single = fread(fidsingle,timestep*freq*nx*ny,'uint16');
    single = reshape(single,nx,ny,timestep*freq);
    for l=1:timestep*freq % iterates every frame within the data block
        D(:,:,l)=(single(:,:,l) - demin_mean(:,:,i))./(ref_mean(:,:,i)*...
            ratio-demin_mean(:,:,i));
        sigma=sign(D(:,:,l));
        D(:,:,l) = D(:,:,l).*(onet+sigma)/2;
        for m=1:nx      %creates the map of mins & maxes for set i
            for n=1:ny
                if mask(m,n)==1
                    %iterates all mesh points for the pdf
                    % is needed on loop prior to the actual PDF loop, same as
                    % D_mean and D_std
                    if D(m,n,l) > D_max(m,n,i);
                        D_max(m,n,i)=D(m,n,l);
                    else
                        if D(m,n,l) < D_min(m,n,i);
                            D_min(m,n,i)=D(m,n,l);
                        end
                    end
                end
            end
        end
    end
    D_meantemp=D_meantemp+sum(D,3);
    save(sprintf('%s\\%s\\D_set%02d_block%02d',pwd,Dfolder,i,j),'D')
end
D_min(:,:,i)=D_min(:,:,i).*mask;
D_max(:,:,i)=D_max(:,:,i).*mask;
D_mean(:,:,i)=D_meantemp.*mask/(freq*singletime);
save(sprintf('%s\\%s\\D_max',pwd,Dfolder),'D_max');
save(sprintf('%s\\%s\\D_min',pwd,Dfolder),'D_min');
save(sprintf('%s\\%s\\D_mean',pwd,Dfolder),'D_mean');

% %% plot
% createfiguremean(D_mean(:,:,i))
% saveas(gcf,sprintf('%s\\%s\\D_mean_%d',pwd,Dfolder,i),'jpg');
% close gcf
%% clean up
%free up memory & keep worspace clean
clear D_meantemp
clear single
clear D 
clear onet
clear sigma
clear D_min
clear D_mean
clear D_max
clear ref_mean
clear demin_mean
clear start stop

