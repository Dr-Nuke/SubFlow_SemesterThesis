%  This scriptfile and its subscripts generate the dimensionless mixing
%  scalar DMS, any derived quantities and plots them in the end.

%% Authors notes
% <taken care of>

%% manual
% 1) place this and all the other m files in the same folder as the data 
% from the SUBFLOW facility (needs to be converted for matlab).

% 2) enter their filenames into names.m. make sure you do this in the
% correct order.

% 3) specify any variable throughout the scripts or leave them as default.

% 4) run this file via F5. It takes several hours to do the calculations.
% results are saved into the corresponding files.

% 5) to look at the files, make use of the skptLoader.m and skptPlotter.m
% scripts for loading ointo workspace and plotting.



%% known bugs & todo list

% i used some times "xxxtime*timestep". i think it should be
% "xxxtime/timestep". i never noticed because timestep=1 since ever

% need to do a proper handling of velx etc: saving, clearing and loading

% vary the timestep and stepsize value and compare results

% fft

%% prerequisites
initialize=1;
if initialize==1    % this if is meant for if you run the code for hours,
                    % get an error, and wantt to continue after debugging
                    % without deleting all your previous results
    clc
    clear all
    close all
    format compact;
    ticstart = tic;
    
    %% data names
    [singlenames,deminnames,refnames]=names(0); 
    % makes data filenames available, caution: bugs are hardly recognised
    % in the running programm, but deliveres wrong data.                                                
    
    %% mask to cut out the fuelrod area (hardcoded ascii file)
    maskimport;
    
    %% settings & constants
    singletime = 20;% time of measurement
    demintime = 5;  % time of calibration measurement
    reftime = 5;    % of the reference measurement
    timestep = 1;   % how much should be read in at the time, for memory 
                    % reasons (in sec)
    nx = 128;       % twice the sensor x resolution (2 sensors next to each 
                    % other)
    ny = 64;        % y resolution
    freq = 2500;    % sampling frequency
    ratio = 787.63; % ratio of the 
    
    %% data import & processing, mean value
    numberofsets=12 ;   % one set consits of one demin, reference and 
                        % single file each
    numberofblocks=singletime/timestep; % by how many do i devide the data
end

%% Part 1: computing the DMS
% (including the min max mean values)
Dfolder='testDMinMaxMean'; %name of the folder to be used for part 1 saves
 
if 0 % use false if the values already exist
    for i = 1:numberofsets % this for loop iterates over the 
                % various sets of inputfiles, see calibration.bat
        if ~exist(Dfolder,'dir')            
            mkdir(Dfolder);
        end
        disp(sprintf('file set %d of %d',i,numberofsets));
        fiddemin=fopen(char(deminnames(i)));    % see help fopen
        fidref=fopen(char(refnames(i)));
        fidsingle=fopen(char(singlenames(i)));
        
        %% demin and reference processing
        skptdeminimport 
        skptrefimport
        %% D_mean, D_min and D_max
        skptDMeanMinMax
    end
end

%% Part 2: computing derived quantities 
STDfolder='testSTD';  % name of the folders to be used for saves
PDFfolder='testPDF'
if 1
	   VELfolder='testvelocitysave'; % name of the folder to be used 
                                          % for saves
           mkdir(VELfolder);
           %NEW LINES ADDED ON 8.2.2012
           load(sprintf('%s\\%s\\D_mean',pwd,Dfolder),'D_mean');
           %
for i = 11:11 %for i= 1:numberofsets  %12 sets of measurements
    %% STD
    %skptStd
    
    %% PDF
    %binnumber=50;       % # of bins
    %skptPDF
    
    %% fft
    %Dpoint =[10,20;    %  meshpoints for the time signal for fft
    %        28,24];    %  enter the coordinates you see in the figures
    %                   %  [x1 y1 x2 y2 x3 y3 ...]
    %skptfft  %requires Dpointtime to be existent

    %% cross correlations
    %     requires a set of data labelled 'D', so by default the 20th data
    %     this is where i need to fix loading & saving etc, for example
    %     package
    %CORRfolder='corr';
    %fp=[17 24 0 24 24 0];  % fastpoints [x1 x2 #x y1 y2 #y],
                            % points (x1:x2,y1:y2) will be referencepoints 
                            % for covariancees
    %skptcrosscorelations

    %% lateral velocity
    if 1
        stepsize=0.5;% fraction of a block of D data that is 
                     % used for each xcov calculation. =>resizes the length
                     % of the time signal that is used for xcov. e.g.
                     % 0.1 means that a block of D is divided in 1/0.1=10
                     % subblocks.
                     % attention: it is believed to dramatically increase
                     % computation time if this number is reduced
        maxlag=125;  % enter the lag in frames.
                     % 15mm/(0.8m/s)=18 ms so i pick 50ms as maxlag,
                     % 0.05s*2500frames/s= 125 frames
        xcovcoeffthreshold=0.5; % sets a threshold below which points 
                                % are omitted
        saveLAT=1;   % see saveD above (skptDMeanMinMax)
        skptLateral  % if condition temporarly changed
    end
end
end
 
%% Part 3: plotting of desired quantities
for i=1:numberofsets



end

    %% Analysis on lateral velcity
%     VELfolder='velocitysave05_tets';
%     saveVELSTD=0;
%     skptLATanalysis

    %% coefficient of variation
    %saveCoefOfVar=1
    %skptcoeffofvariation
      
     %% reduced velocities
     
     %% 3d velocity plots
     
     %% 

%% 
disp(sprintf(' %d seconds total', toc(ticstart)))