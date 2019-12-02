function [Dmean,Dstd] = wmsanalysis(dname,deminname,tracername,ratio,time,mask)
% Example calling:
%[Dmean,Dstd] = wmsanalysis('single_1_1','demin_1_1_1','reference_1_2',787.63,20,mask);
%close all
caltime = 5;
freq=2500;
nx = 128;
ny = 64;
tic;
filetype = '.dat';
coef=zeros(nx,ny);

%% Reading of the mask
%fid_mix = fopen('mask.asc','r');
%mask = fscanf(fid_mix,'%g',[nx ny]);
for i=1:nx
for j=1:ny
    if mask(i,j)==0
    mask(i,j)=NaN;
    end
end
end
%fclose(fid_mix);

%% Calibration file for Demineralized water
filename = strcat(deminname, filetype);
demin_mean = zeros(nx,ny);
fid = fopen(filename);

for kk=1:caltime
demin = fread(fid,freq*nx*ny,'uint16');
demin = reshape(demin,nx,ny,freq);
for k=1:freq
        demin_mean = demin_mean+demin(:,:,k);
end
end
fclose(fid);
demin_mean = demin_mean/(caltime*freq);

%% Calibration file for Tracer
filename = strcat(tracername, filetype);
reference_mean = zeros(nx,ny);
fid = fopen(filename);

for kk=1:caltime
reference = fread(fid,freq*nx*ny,'uint16');
reference = reshape(reference,nx,ny,freq);
for k=1:freq
        reference_mean = reference_mean+reference(:,:,k);
end
end
fclose(fid);
reference_mean = reference_mean/(caltime*freq);

clearvars reference demin
%% Reading of the measurement file
filename = strcat(dname, filetype);
Dmean = zeros(nx,ny);
Dstd = zeros(nx,ny);
fid = fopen(filename);

for kk=1:time
D = fread(fid,freq*nx*ny,'uint16');
D = reshape(D,nx,ny,freq);

%% Calculation of Dimensionless mixing scalars
signmatrix = zeros(nx,ny);
onet = ones(nx,ny);

for k=1:freq;
        D(:,:,k) = (D(:,:,k)-demin_mean)./(reference_mean*ratio-demin_mean);
        signmatrix = sign(D(:,:,k));
        D(:,:,k) = D(:,:,k).*((onet+signmatrix)/2);
        Dmean = Dmean+D(:,:,k);
end
end
fclose(fid);
Dmean = (Dmean/(time*freq)).*mask;

%%Calculation of standard deviation
fid = fopen(filename);

for kk=1:time
D = fread(fid,freq*nx*ny,'uint16');
D = reshape(D,nx,ny,freq);

signmatrix = zeros(nx,ny);
onet = ones(nx,ny);

for k=1:freq;
        D(:,:,k) = (D(:,:,k)-demin_mean)./(reference_mean*ratio-demin_mean);
        signmatrix = sign(D(:,:,k));
        D(:,:,k) = D(:,:,k).*((onet+signmatrix)/2);
        Dstd = Dstd+(D(:,:,k)-Dmean(:,:)).^2;
end
end
fclose(fid);

Dstd = sqrt((Dstd/(time*freq-1)).*mask);

%Coefficient of variation
%Co=std/mean;




% %% Standard deviation of Dimensionless mixing scalar
% Dsumsqr = zeros(nx,ny);
% for k=1:time*freq
%     Dsumsqr=Dsumsqr+(D1(:,:,k)-Dmean).*(D1(:,:,k)-Dmean);
% end
% Dstd = sqrt(Dsumsqr/(time*freq));
% 
% %% Calculation of cross-correlation map
% if (crossi==1)
% 
% % Setting of the fixed point
% ref = reshape(D1(17,24,:),1,time*freq);
% 
% for i=1:nx
%  for j=1:ny
%         loc=reshape(D1(i,j,:),1,time*freq);
%         coef(i,j)=fastcorr(ref,loc);
%  end
% end
% coef = fliplr(flipud(coef));
% 
% h1 = figure('Position', [100 100 800 800]);
% subplot('Position',[0.05 0.05 0.89 0.89]);
% [c,h] = contourf(coef,30,'LineStyle',':','LineWidth',1);
% colorbar
% axis square
% set(gca, 'fontsize',26);
% axis off
% colormap(hot)
% end
% 
% %% PDFs from dimensionless mixing scalars
% figure
% hold on
% 
% MarkerEdgeColors=['r','b','g','c','m','k','m','c','g','b','r'];
% linestyles = cellstr(char('-','-','-','-','-','-','-.','-.','-.','-.','-.'));
% j=1;
% 
% for i=4:45
% sig=D1(24,i,10:length(D1));
% sig=reshape(sig,1,length(sig));
% [f,xi]=ksdensity(sig);
% summa=sum(f);
% f1=f/summa;
% if (i==4)|(i==8)|(i==12)|(i==16)|(i==20)|(i==24)|(i==29)|(i==33)|(i==37)|(i==41)|(i==45)
% plot(xi,f1, [linestyles{j} MarkerEdgeColors(j)],'LineWidth',2)
% j=j+1;
% end
% end
% legend('Point 11','Point 21','Point 31','Point 41','Point 51','Injection point','Point 52','Point 42','Point 32','Point 22','Point 12')
% set(gca, 'fontsize',16);
% hold off
% 
% %% Creation of movie frames
% if (vidset==1)
% % Number of frames
% N = 500;
% T = 0;
% T0 = 10/freq;
% frame = zeros(64,64,N);
% h1 = figure('Position', [100 100 800 800]);
% subplot('Position',[0.05 0.05 0.95 0.8]);
% %subplot('Position',[left bottom width height])
% set(gca,'NextPlot','replacechildren')
% 
% annotation(h1,'textbox',[0.15 0.87375 0.30 0.1],'Interpreter','none',...
%     'String',{'SUBFLOW',dname},...
%     'FontSize',25,...
%     'FontName','Arial',...
%     'FitBoxToText','off',...
%     'LineStyle','none');
% 
% colorbar('FontSize',20);
% colormap(hot)
% set(gca, 'fontsize',20);
% axis square
% 
% for k=10:(10+N)
% frame(:,:,k-9) = D1(:,:,k);
% frame(:,:,k-9) = reshape(rot90(frame(:,:,k-9),2),nx,ny).*mask;
% end
% maxi = max(max(max(frame)));
% caxis([0 maxi])
% 
% for k=10:(10+N)
% [c,h] = contourf(frame(:,:,k-9),30);
% set(h,'edgecolor','none')
% 
% % Create colorbar
% T = T0+k/freq;
% title(['Time ' num2str(sprintf( ' %10.4f ' , T )) 's'],'FontSize',25,'FontName','Arial');
% 
% saveas(h1,[pathout,dname,'_' num2str(k)],'png');
% end
% close all
% end

%%

toc;

end
