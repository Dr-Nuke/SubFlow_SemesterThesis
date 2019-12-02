%% initialisation:

subblocks=singletime/(timestep*stepsize); % number of the blocks iterated with n & o loop below
xcovcoeff=zeros(1,2*maxlag-1);          % temporary storage
lags=zeros(1,2*maxlag-1);               % also temporarystorage
xcovcoeff2=zeros(nx/2,ny,subblocks);    % will contain the max coeff of
                                        % vincinity of (j,k) for all
                                        % subblocks
xcovlag=zeros(nx/2,ny,subblocks);       % partner of xcovcoeff2
xcovcoordx=NaN(nx/2,ny,subblocks);      % not introduced as 0 because 0 is a reasonable value
xcovcoordy=NaN(nx/2,ny,subblocks);
D=zeros(nx,ny,timestep*freq);           % data
velocityx=zeros(nx/2,ny,subblocks);     % the various velocity profiles
velocityy=zeros(nx/2,ny,subblocks);
velocityz=zeros(nx/2,ny,subblocks);
velocitylat=zeros(nx/2,ny,subblocks);
velocitytot=zeros(nx/2,ny,subblocks);
velocityzcounter=zeros(nx/2,ny);        % needed for the mean calculation

%% variables that see all i
if  exist(sprintf('%s\\%s\\veltot.mat',pwd,VELfolder),'file')  % these variables are for all i of one set of data
    load(sprintf('%s\\%s\\velx',pwd,VELfolder),'velx');
    load(sprintf('%s\\%s\\vely',pwd,VELfolder),'vely');
    load(sprintf('%s\\%s\\velz',pwd,VELfolder),'velz');
    load(sprintf('%s\\%s\\vellat',pwd,VELfolder),'vellat');
    load(sprintf('%s\\%s\\veltot',pwd,VELfolder),'veltot');
else
    velx=zeros(nx/2,ny,numberofsets);
    vely=zeros(nx/2,ny,numberofsets);
    velz=zeros(nx/2,ny,numberofsets);
    vellat=zeros(nx/2,ny,numberofsets);
    veltot=zeros(nx/2,ny,numberofsets);
end
%flowvel=0.8;    % the average velocity of the flow single_1_1
%flowvel=1.2;    % not needed in calculations but there to check results
j=1;
k=1;
l=1;
m=1;
NN=0; %Neigbourhood used in the calculations
distz=0.015;     % distance of the 2 sensors
mesh=0.002125;  % resolution of the sensor
%% calculation
if saveLAT==1  % either compute or load the results
    for n=1:numberofblocks  % corresponds to the number of files for D
        load(sprintf('%s\\%s\\D_set%02d_block%02d',pwd,Dfolder,i,n),'D'); %load data
        for o=1:timestep/stepsize   % subblocks of one D block of data
            start = round(((o-1)*stepsize*freq)+1); %subblock boundaries in terms of frames
            stop = round(o*stepsize*freq);
            disp([start stop]);
            nn= (n-1)*timestep/stepsize+o;  % continous counter. counts all subblocks for one i
            
            %% cross correlation coefficients and lags
            disp('lateral');
            for j=1:nx/2       % iterates all pixel in x and
                disp(sprintf('n=%d o=%d nn=%d j=%d k=%d l=%d m=%d i=%d %d %d %d %d ',n,o,nn,j,k,...
                    l,m,i,size(mask,2),mask(1,1),nx,ny)) % debug tool
               for k=1:ny     % in y direction
                    if (mask(j,k)==1) && (D_mean(j,k,i)>0.00001)  % no need to check inside the fuelrods...
                        for l = max(1,j-NN):min(nx/2,j+NN) % goes from j-NN to j+NN except at
                                                         % the boundaries of the array
                            for m = max(1,k-NN):min(ny,k+NN) %  l for x, m for y
                                cntrl = mask(l,m);
                                if cntrl==1;  %no need to check inside the fuelrods...
                                     
                                        [xcovcoeff,lags]=xcov(squeeze(D(j+nx/2,k,start:stop)),...
                                                              squeeze(D(j,k,start:stop)),maxlag,'coeff');
                                        if max(xcovcoeff) > xcovcoeffthreshold
                                        
                                        T = lags./freq;   % Assuming cross covariance coefficients as PDF of flying time(?). 
                                                          % Then PDF of velocity can be derived according to relationship v = distz/time
                                                          % ATTENTION: average velocity obtained in this way is about 0.45m/s... wrong...                                                          
                                      
                                       % flow direction is known, it can be presumed lags coresponding to max xcovcoeff are always positive
                                        for p = 1: 2*maxlag+1  
                                            if T(p) <= 0
                                                T(p) = NaN;
                                            end
                                        end
                                        pdfT = xcovcoeff;
                                        v = distz./T;
                                        pdfv = pdfT.*T'.*T'*distz;
                                        
                                        vmost = v(find(pdfv == max(pdfv)));  % most probable velocity according to PDF of velocity
                                        velocitytot(j,k,nn) = vmost;
                                        end
%                                         figure
%                                         plot(T,pdfT)
%                                         figure
%                                         plot(v,pdfv)

                                                          
%                                         if max(abs(xcovcoeff))> xcovcoeff2(j,k,nn)
%                                             xcovcoeff2(j,k,nn)=max(abs(xcovcoeff));
%                                             xcovlag(j,k,nn)=find(abs(xcovcoeff) == max(abs(xcovcoeff)))-(maxlag+1);
%                                             xcovcoordx(j,k,nn)=l-j;
%                                             xcovcoordy(j,k,nn)=m-k;
%                                         end
                                         
                                    
                                end
                            end
                        end
                    end
                end
            end
            %now, regarding the vincinity of each point in the first sensor, we
            %have maps of the maximum crosscorrelation and the corresponding
            %lags and lateral shifts
            
            %% creating the velocities
%             disp('z-velocity');
%             for j=1:nx/2       % iterates all pixel in x and
%                 for k=1:ny     % in y direction
%                     
%                     if mask(j,k)==1        % no need to check inside the fuelrods...
%                         if (xcovcoeff2(j,k,nn)>xcovcoeffthreshold)  && (xcovlag(j,k,nn)~=0)
%                                         velocityz(j,k,nn)=distz/(xcovlag(j,k,nn)/freq);
%                                         velocityzcounter(j,k)=velocityzcounter(j,k)+1;
%                                         velocitytot(j,k,nn)=sqrt(distz^2+(xcovcoordx(j,k,nn)*mesh)^2 ...
%                                         +(xcovcoordy(j,k,nn)*mesh)^2)/(xcovlag(j,k,nn)/freq);
%                                         velocityx(j,k,nn)=mesh*xcovcoordx(j,k,nn)/(xcovlag(j,k,nn)/freq);
%                                         velocityy(j,k,nn)=mesh*xcovcoordy(j,k,nn)/(xcovlag(j,k,nn)/freq);
%                                         velocitylat(j,k,nn)=abs(sqrt(velocityy(j,k,nn)^2+velocityx(j,k,nn)^2));
%                         end
%                     end
%                 
%                 end
%             end                   
%            
       

            %% print & save intermediate results 
            % may produce rediculus anounts of files for small "stepsize"
%             %% printing and saving
%             %if nn==40
%             createfigurevel(velocityx(:,:,nn),1);
%             saveas(gcf,sprintf('%s\\%s\\velocityx_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%             close gcf
%             %title('velocityx')
%             createfigurevel(velocityy(:,:,nn),2);
%             saveas(gcf,sprintf('%s\\%s\\velocityy_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%             close gcf
%             %title('velocityy')
%             createfigurevel(velocityz(:,:,nn),3);
%             saveas(gcf,sprintf('%s\\%s\\velocityz_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%             close gcf
%             %title('velocityz')
%             createfigurevel(velocitylat(:,:,nn),4);
%             saveas(gcf,sprintf('%s\\%s\\velocitylat_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%             close gcf
%             %title('velocitylat')
%             createfigurevel(velocitytot(:,:,nn),5);
%             saveas(gcf,sprintf('%s\\%s\\velocitytot_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%             close gcf
%             %title('velocitytot')
%             %figure
%             % end
        end
    end
        
    %% the velocity values: creation, saving, plotting
    
%     velx(:,:,i)=sum(velocityx,3)./velocityzcounter;     %covers all i  % =mean(velocityx,3)
%     save(sprintf('%s\\%s\\velx',pwd,VELfolder),'velx');
%     save(sprintf('%s\\%s\\velocityx_%d',pwd,VELfolder,i),'velocityx');
%         
%     vely(:,:,i)=sum(velocityy,3)./velocityzcounter;
%     save(sprintf('%s\\%s\\vely',pwd,VELfolder),'vely');
%     save(sprintf('%s\\%s\\velocityy_%d',pwd,VELfolder,i),'velocityy');
%     
%     velz(:,:,i)=sum(velocityz,3)./velocityzcounter;
%     save(sprintf('%s\\%s\\velz',pwd,VELfolder),'velz');
%     save(sprintf('%s\\%s\\velocityz_%d',pwd,VELfolder,i),'velocityz');
%     
%     vellat(:,:,i)=sum(velocitylat,3)./velocityzcounter;
%     save(sprintf('%s\\%s\\vellat',pwd,VELfolder),'vellat');
%     save(sprintf('%s\\%s\\velocitylat_%d',pwd,VELfolder,i),'velocitylat');
    
    veltot(:,:,i)=sum(velocitytot,3)./velocityzcounter;
    save(sprintf('%s\\%s\\veltot',pwd,VELfolder),'veltot');
    save(sprintf('%s\\%s\\velocitytot_%d',pwd,VELfolder,i),'velocitytot');
    
%     createfigurevel(velx(:,:,i),1);
%     saveas(gcf,sprintf('%s\\%s\\velocityxfin_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%     close gcf
%     createfigurevel(vely(:,:,i),2);
%     saveas(gcf,sprintf('%s\\%s\\velocityyfin_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%     close gcf
%     createfigurevel(velz(:,:,i),3);
%     saveas(gcf,sprintf('%s\\%s\\velocityzfin_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%     close gcf
%     createfigurevel(vellat(:,:,i),4);
%     saveas(gcf,sprintf('%s\\%s\\velocitylatfin_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
%     close gcf
    createfigurevel(veltot(:,:,i),5);
    saveas(gcf,sprintf('%s\\%s\\velocitytotfin_%d_%d_%d',pwd,VELfolder,i,nn,nn),'jpg');
    close gcf
    
    save(sprintf('%s\\%s\\velocityzcounter_%d',pwd,VELfolder,i),'velocityzcounter');
    save(sprintf('%s\\%s\\xcovcoordx_%d',pwd,VELfolder,i),'xcovcoordx');
    save(sprintf('%s\\%s\\xcovcoordy_%d',pwd,VELfolder,i),'xcovcoordy');
    save(sprintf('%s\\%s\\xcovlag_%d',pwd,VELfolder,i),'xcovlag');
    save(sprintf('%s\\%s\\xcovcoeff_%d',pwd,VELfolder,i),'xcovcoeff2');
end %end of the saveLAT==1 loop
if saveLAT==0
    load(sprintf('%s\\%s\\velx',pwd,VELfolder),'velx');
    load(sprintf('%s\\%s\\vely',pwd,VELfolder),'vely');
    load(sprintf('%s\\%s\\velz',pwd,VELfolder),'velz');
    load(sprintf('%s\\%s\\vellat',pwd,VELfolder),'vellat');
    load(sprintf('%s\\%s\\veltot',pwd,VELfolder),'veltot');
    % maybe the loading should be done in the script that needs this
end
clear D
clear velocityx
clear velocityy
clear velocityz
clear velocitylat
clear velocitytot
clear xcovcoeff
clear lags
clear xcovcoeff2