%% initialisation
if ~exist(PDFfolder,'dir')
    mkdir(PDFfolder);
end
if ~exist(sprintf('%s\\%s\\pdf_%02d.mat',pwd,PDFfolder,i),'file')
    pdf=zeros(numberofsets,nx,ny,binnumber);     % will contain PDF data
    pdfx=zeros(numberofsets,nx,ny,binnumber);    % pdf x-axis
                            
else
    load(sprintf('%s\\%s\\pdf_%02d',pwd,PDFfolder,i),'pdf')
    load(sprintf('%s\\%s\\pdfx_%02d',pwd,PDFfolder,i),'pdfx')    
end
pdfxbackup=zeros(1,nx,ny,2);    % is needed for an +-inf issue, see below
pdftemp=zeros(1,binnumber);     % will contain PDF data  of one loop
D=zeros(nx,ny,timestep*freq);   % memory allocation
disp(sprintf('creating PDF set %d x axes',i));
load(sprintf('%s\\%s\\D_min',pwd,Dfolder),'D_min');
load(sprintf('%s\\%s\\D_max',pwd,Dfolder),'D_max');,
%% PDF Xaxis

for m=1:nx %  these two forloops create the x axis' for all PDFs
    for n=1:ny
        cntrl = mask(m,n);  % no need to calculate in the fuel rods
        if cntrl == 1
            pdfx(i,m,n,1:binnumber)=linspace(D_min(m,n,i),D_max(m,n,i),...
                binnumber); % the bin-boundaries for the pdf
        end
    end
end

pdfxbackup(1,:,:,1)=pdfx(i,:,:,1);
pdfxbackup(1,:,:,2)=pdfx(i,:,:,end);
                % the plot-x-axis cant handle the inf values 
                % but i want them for the pdf. afterwards i reset the
                % values
pdfx(i,:,:,1)=-inf*ones(nx,ny);
pdfx(i,:,:,end)=inf*ones(nx,ny);    % to include ALL range 
                                    % maybe obsolete due to D_minmax?

%% PDF
disp(sprintf('skpPDF: processing D frame _ to _ out of %d',freq*singletime));
status = fseek(fidsingle,0,'bof'); %set the pointer to beginning of file
for j = 1:singletime/timestep   % iterates the blocks of data 
    start = ((j-1)*timestep*freq)+1;
    stop = j*timestep*freq;
    disp([start stop]);
    load(sprintf('%s\\%s\\D_set%02d_block%02d',pwd,Dfolder,i,j),'D')

    for m=1:nx             % m is the x coordinate 
        for n=1:ny
            cntrl = mask(m,n);
            if cntrl == 1   % indicates wether the point is outside 
                            % the fuelrods. saves some time
                pdftemp=histc(D(m,n,:),pdfx(i,m,n,:));
                pdf(i,m,n,:)=reshape(pdftemp,1,1,1,binnumber)+pdf(i,m,n,:); % NOT normalizing here
            end
        end
    end
end

pdfx(i,:,:,1)=pdfxbackup(1,:,:,1);
pdfx(i,:,:,end)=pdfxbackup(1,:,:,2);
save(sprintf('%s\\%s\\pdf_%02d',pwd,PDFfolder,i),'pdf')
save(sprintf('%s\\%s\\pdfx_%02d',pwd,PDFfolder,i),'pdfx')


% pdf=pdf/(freq*singletime); %i implemented this at the very end of
% developping time but im convinced its wrong. now (without it) the sum
% of one pdf is50000 %[2nd edit: that seems to be wrong too....]

%% cleanup
clear pdf
clear pdfx
clear pdfxplot
clear pdftemp
clear start stop
clear pdfbackup
clear D_min D_max
clear pdfxbackup
clear D


% %% plot the PDFs
%   
% cc=jet(128);
% figure 
% hold on
% linecolor=hot(11);
% m=0; %graph counter
% for j=64:128 
%     for k=9
%        
%        
%         if 1 ;     % indicates wether the point is outside the fuelrods
% 
%             plot(reshape(pdfx(i,j,k,:),1,binnumber),reshape(pdf(i,j,k,:),1,binnumber),'color',cc(j,:));
%             m=m+1;
%            % xlim([0 0.002])
%             title(j)
%         end
% 
%     end
% end
