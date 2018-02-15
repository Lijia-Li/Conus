function runAnalysis3(handles)
%%%%% VERSION 3.1 January 2014
%%%%% For Windows/Mac/Unix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get input parameters from GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FilesPath: The file to analyze
filesPath = get(handles.FilesPath,'String');
% outputPath: Where output information will be written
outputPath = get(handles.OutputPath,'String');
% fishSet: The fish to be analyzed (see processing section)
outputFile = fullfile(outputPath,get(handles.FileName,'String'));
ppoTimestr = get(handles.ppoTime,'String');
ppoTime = str2num(get(handles.ppoTime,'String'));
noiseThresh = str2double(get(handles.noiseThresh,'String'));
plotGMVOT = get(handles.plotGMVOT,'Value');
plotGPTOT = get(handles.plotGPTOT,'Value');
plotNoise = get(handles.plotNoise,'Value');
plotHeatmap = get(handles.plotHeatmap,'Value');
plotIntensities = get(handles.plotIntensities,'Value');
plotHistograms = get(handles.plotHistograms,'Value');
plotPPO = get(handles.plotPPO,'Value');
plotColormap = get(handles.plotColormap,'Value');
plotGMQOT = get(handles.plotGMQOT,'Value');
plotGQPTOT = get(handles.plotGQPTOT,'Value');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Input and Output Information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numMovies = 0;
if (isempty(filesPath) || isempty(outputPath) || strcmp(filesPath,'No files selected...') || strcmp(outputPath,'No directory selected...'))
    errordlg('You did not specify input file(s) or an output directory');
    set(handles.goBut,'Enable','on');
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load mat file and further process parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(filesPath);

% strcmp compares two strings.  If the files to analyze field contains
% "[:]" there is only one analysis group, 
% otherwise, this part of the routine takes the input from the interface
% and uses it to generate well lists that belong to each experiemntal group
% for example if two groups of 48 wells are chosen fishSet = [1x48 double] [1x48 double] 
% fishSet is a cell format
% in example above, fish set contains two cells, each with 1, 2, 3, 4 ...
% and 49, 50, 51
% which specifies the well numbers that belong to each group

if(strcmp(get(handles.fishSet,'String'),'[:]'))
    fishSet{1} = [1:size(fishDistances,1)];
else 
    fishSet = str2cell(get(handles.fishSet,'String'));
end

allfish = cell2mat(fishSet);

if(strcmp(get(handles.ppoTime,'String'),'[:]'))
    ppoTime = [1:size(fishDistances,2)];
end

fprintf('\nAnalyzing data and generating figures, please wait...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set parameters for detection and display
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smoothFactor = 50;   % the data calculated below will be applied with a smooth filter
emptyWellThresh = 5; % min % of frames with NOEP before well is considered empty 
maxNoiseThresh = noiseThresh*100; % max % of frames that TMOEP or NOEP can be detected before well is thrown out


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Convert from pixels/frame to mm/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('frameRate')
    frameRate = 2;
end
[numWells,numFrames] = size(fishDistances);
%To add functionality for other wells, put the plate well number in the
%first set, and the corresponding well diameter in the same spot of the second 
wellDiameterConv = containers.Map({96,48,40,24,12,6},{6.78,10.5,9.8,15.62,22.1,34.8});  %% 7/16 are optionally used for 96/24 wells
%pix/frame * diameter(mm)/2*radius(pix) * frameRate frames/second
if ~exist('mmConv')
    disp(wellDiameterConv(numWells))
    mmConv = (wellDiameterConv(numWells)/(2*unscaledRadius))*frameRate; 
end
allfishVelocities = fishDistances*mmConv; %mm/s

allfishQuants = fishQuants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare output structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GroupOut = {'Wells','n', 'Mean Velocity (mm/s)', 'Mean Velocity SE', ...
    'Active Velocity (mm/s)', 'Active Velocity SE', '% Time Moving',...
    '% Time Moving SE', 'Active Duration (s)', 'Active Duration SE', ...
    'Rest Duration (s)', 'Rest Duration SE'};
IndividOut = {'Set','Well','Mean Velocity (mm/s)', 'Mean Velocity SD','Active Velocity (mm/s)','% Time Moving','Active Duration (s)', 'Rest Duration(s)'};
GroupOutQ = {'Wells','n', 'Mean Quant (/s)', 'Mean Quant SE', ...
    'Active Quant (/s)', 'Active Quant SE', '% Time Moving',...
    '% Time Moving SE', 'Active Duration (s)', 'Active Duration SE', ...
    'Rest Duration (s)', 'Rest Duration SE'};
IndividOutQ = {'Set','Well','Mean Quant (/s)', 'Mean Quant SD','Active Quant (/s)','% Time Moving','Active Duration (s)', 'Rest Duration(s)'};

% To properly process 
firstTime = 1;
for setNum = 1:length(fishSet)
    %Get the current fish group
    fish = fishSet{setNum};
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Assess well usability
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate a list of all possible usable wells
    list = 1:numWells;
    % Find wells with acceptable noObjectError percent (NOEP) and
    % tooManyObjectError percent (TMOEP)
    NOEP = noObjectError./numFrames*100;
    okNOEP = intersect(fish,list(NOEP < maxNoiseThresh));
    TMOEP = tooManyObjectError./numFrames*100;
    okTMOEP = intersect(fish,list(TMOEP < maxNoiseThresh));
    % Find empty wells (those with NOEP > empty threshold) 
    emptyWells = intersect(fish,list(NOEP > emptyWellThresh));
    % Find clean wells (those that have ok NOEP and TMOEP error rates)
    cleanWells = intersect(okNOEP,okTMOEP);
    % Find dirty wells (those that are not clean (high NOEP and TMOEP error rates))
    dirtyWells = setdiff(fish, cleanWells);
    % Find Usable Wells (those that are clean but not empty)
    nonemptyWells = setdiff(fish,emptyWells);
    usableWells = setdiff(cleanWells, emptyWells);
    fishSet{setNum} = usableWells;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Prepare data for analysis 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some functions don't like scalar data (eg. ribbon plots)
    %  A trick to keep singular samples multidimentional
    %  is to simply double the sample to create a group.
    %  This has no effect on mean, std, or other measures, and will show
    %  an interfish variability of 0 (which is accurate).
    %  The only odd effect of this is that the ouput lists the sample twice and
    %  calls the group size 2. The alternative is that  everything below would 
    %  have to be rewritten for the degenerative matrix case (a vector). 
    if(length(usableWells) == 1)
        usableWells(end+1)=usableWells;
    end

    %If no wells are clean, skip this group
    if (isempty(usableWells))
        fprintf('\nWell(s) %s is(are) unusable (>5%% noise, >+/-2SE, empty, or n = 1',num2str(usableWells));
        continue;
    end
    % Save old fishVelocities (for debugging purposes, and to restore when finished)
    oldFishVelocities = allfishVelocities;
    % Store Distances only for usablewells
    fishVelocities = allfishVelocities(usableWells,:);
    [numWells,numFrames] = size(fishVelocities);
    fishQuants = allfishQuants(usableWells,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate Mean & Active Velocities, % Time 
    %%% Active and Group mean wrt Time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (plotGMVOT)
        % GMVOT = group mean velocity over time
        GMVOT = mean(fishVelocities);
        GMVOTsmooth = gaussSmooth(GMVOT,floor(length(GMVOT)/smoothFactor)+1);
    end
    
    if (plotGPTOT)
        % GPTOT = group percent time moving over time according to tracking
        GPTOT = logical(fishVelocities);
        GPTOT = mean(GPTOT);
        GPTOTsmooth = gaussSmooth(GPTOT,floor(length(GPTOT)/smoothFactor)+1);
    end
    
    if (plotGMQOT)
        % GMQOT = group mean quantification over time
        GMQOT = mean(fishQuants);
        GMQOTsmooth = gaussSmooth(GMQOT,floor(length(GMQOT)/smoothFactor)+1);
    end
    
    if(plotGQPTOT)
        % GQPTOT = group percent time moving over time according to quantification
        GQPTOT = logical(fishQuants);
        GQPTOT = mean(GQPTOT);
        GQPTOTsmooth = gaussSmooth(GQPTOT,floor(length(GQPTOT)/smoothFactor)+1);
    end
    
    individMVs = mean(fishVelocities');
    individMQs = mean(fishQuants');    
    individMVSTDs = std(fishVelocities');
    individMQSTDs = std(fishQuants');
    
    %%%% The sum of all velocities divided by the number of velocites > 0 is
    %%%% how active velocity is defined
    individAVs = sum(fishVelocities')./sum(fishVelocities'>0);
    individAQs = sum(fishQuants')./sum(fishQuants'>0);
    %%%% The percent time movement is the number of velocites > 0 over the
    %%%% total number of velocities
    individTPs = (sum(fishVelocities'>0)./numFrames).*100;
    individQTPs = (sum(fishQuants'>0)./numFrames).*100;
    %%%% NOTE: You may want to change the 0's to something else for AV and TP
    %%%% depending on how you define movement and account for noise
    %Remove any NaN or Inf values from active velocity by setting them to 0
    individAVs((or(isinf(individAVs),isnan(individAVs))))=0;
    individAQs((or(isinf(individAQs),isnan(individAQs))))=0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculate burst and rest durations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    individADs(numWells) = 0;
    individRDs(numWells) = 0;
    for x = 1:numWells
        active = regionprops(im2bw(fishVelocities(x,:)),'Area');
        rest = regionprops(not(fishVelocities(x,:)),'Area');
        individADs(x) = mean([active.Area])/frameRate;
        individRDs(x) = mean([rest.Area])/frameRate;
    end   
    individQADs(numWells) = 0;
    individQRDs(numWells) = 0;
    for x = 1:numWells
        Qactive = regionprops(im2bw(fishQuants(x,:)),'Area');
        Qrest = regionprops(not(fishQuants(x,:)),'Area');
        individQADs(x) = mean([Qactive.Area])/frameRate;
        individQRDs(x) = mean([Qrest.Area])/frameRate;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Display group information
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nVideo Frame Rate (frames/sec): %.2f\tConversion Factor: %.2f',frameRate, mmConv);
    fprintf('\n\n_________________ Group %s _________________________________________',num2str(setNum));
    fprintf('\n================ Well Information =================================');
    fprintf('\nEmpty wells:      %s', num2str(emptyWells));
    fprintf('\nDiscarded wells: %s', num2str(setdiff(dirtyWells,emptyWells)));
    fprintf('\nAnalyzed wells:     %s', num2str(usableWells));
    fprintf('\n# of wells analyzed: %s', num2str(length(usableWells)));
    fprintf('\n================ Performance Analysis =============================');
    fprintf('\n      Type                 \tAverage                               ');
    fprintf('\n"No object" errors:        \t%.2f %%', nanmean(NOEP(okNOEP)));
    fprintf('\n"Too many objects" errors: \t%.2f %%', nanmean(TMOEP(okTMOEP)));
    fprintf('\n================ Fish Tracking Activity ===========================');
    fprintf('\n      Type           \tAverage     [Standard Error]   ');
    fprintf('\nMean Velocity:       \t%.2f mm/s\t  [%.2f]', nanmean(individMVs), ste(individMVs));
    fprintf('\nActive Velocity:     \t%.2f mm/s\t  [%.2f]', nanmean(individAVs), ste(individAVs));
    fprintf('\nPercent Time Moving: \t%.1f %%   \t  [%.2f]\n', nanmean(individTPs), ste(individTPs));
    fprintf('\nActive Duration:     \t%.2f sec \t  [%.2f]',nanmean(individADs), ste(individADs));
    fprintf('\nRest Duration:       \t%.2f sec \t  [%.2f]',nanmean(individRDs), ste(individRDs));
    fprintf('\n___________________________________________________________________\n');
    fprintf('\n================ Fish Quantification Activity =====================');
    fprintf('\n      Type           \tAverage     [Standard Error]   ');
    fprintf('\nMean Quant:          \t%.2f /s\t  [%.2f]', nanmean(individMQs), ste(individMQs));
    fprintf('\nActive Quant:        \t%.2f /s\t  [%.2f]', nanmean(individAQs), ste(individAQs));
    fprintf('\nPercent Time Moving: \t%.1f %%   \t  [%.2f]\n', nanmean(individQTPs), ste(individQTPs));
    fprintf('\nActive Duration:     \t%.2f sec \t  [%.2f]',nanmean(individQADs), ste(individQADs));
    fprintf('\nRest Duration:       \t%.2f sec \t  [%.2f]',nanmean(individQRDs), ste(individQRDs));
    fprintf('\n___________________________________________________________________\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Store info for file output and graphing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GroupOut = cat(1,GroupOut,{usableWells,length(usableWells),mean(individMVs),ste(individMVs),...
        mean(individAVs), ste(individAVs), mean(individTPs), ste(individTPs),...
        mean(individADs), ste(individADs), mean(individRDs), ste(individRDs)});
    GroupOutQ = cat(1,GroupOutQ,{usableWells,length(usableWells),mean(individMQs),ste(individMQs),...
        mean(individAQs), ste(individAQs), mean(individQTPs), ste(individQTPs),...
        mean(individQADs), ste(individQADs), mean(individQRDs), ste(individQRDs)});
    if (firstTime)
        if (plotGMVOT)
            GMVOTplot = GMVOT;
            GMVOTsmoothplot = GMVOTsmooth;
        end
        if (plotGPTOT)
            GPTOTplot = GPTOT;
            GPTOTsmoothplot = GPTOTsmooth;
        end
        if (plotGMQOT)
            GMQOTplot = GMQOT;
            GMQOTsmoothplot = GMQOTsmooth;
        end
        if (plotGQPTOT)
            GQPTOTplot = GQPTOT;
            GQPTOTsmoothplot = GQPTOTsmooth;
        end
        usedWells = usableWells;
        usedGroups = strcat(sprintf('Group %2d', setNum));
    else
        if (plotGMVOT)
            GMVOTplot = cat(1,GMVOTplot,GMVOT);
            GMVOTsmoothplot = cat(1,GMVOTsmoothplot,GMVOTsmooth);
        end
        if (plotGPTOT)
            GPTOTplot = cat(1,GPTOTplot,GPTOT);
            GPTOTsmoothplot = cat(1,GPTOTsmoothplot,GPTOTsmooth); %calculate the smoothed data for overlay display
        end
        if (plotGMQOT)
            GMQOTplot = cat(1,GMQOTplot,GMQOT);
            GMQOTsmoothplot = cat(1,GMQOTsmoothplot,GMQOTsmooth);
        end
        if (plotGQPTOT)
            GQPTOTplot = cat(1,GQPTOTplot,GQPTOT);
            GQPTOTsmoothplot = cat(1,GQPTOTsmoothplot,GQPTOTsmooth); %calculate the smoothed data for overlay display
        end 
        usedWells = cat(2,usedWells,usableWells);
        usedGroups = cat(1,usedGroups,strcat(sprintf('Group %2d', setNum)));
    end
    for i = 1:length(usableWells)
        IndividOut = cat(1,IndividOut,{setNum,usableWells(i),individMVs(i),individMVSTDs(i),individAVs(i),individTPs(i),individADs(i),individRDs(i)});
        IndividOutQ = cat(1,IndividOutQ,{setNum,usableWells(i),individMQs(i),individMQSTDs(i),individAQs(i),individQTPs(i),individQADs(i),individQRDs(i)});
    end

    %errorPlot = NOEP(usedWells)+TMOEP(usedWells);
    %errorPlot = NOEP(:)+TMOEP(:);

    fishVelocities = oldFishVelocities;
    [numWells,numFrames] = size(fishVelocities);
    firstTime = 0;

    % clean some values
    clear individMVs individAVs individTPs individADs individRDs usableWells;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Write output files
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Write the group excel file mean and standard errors
% GroupOut = cellfun(@num2str,GroupOut,'UniformOutput',false);
% GroupOutQ = cellfun(@num2str,GroupOutQ,'UniformOutput',false);
% %%% Write the individual excel file mean and stds
% temp = cellfun(@num2str,IndividOut,'UniformOutput',false);
% tempQ = cellfun(@num2str,IndividOutQ,'UniformOutput',false);
% 
% xlsfilename = strcat(outputFile(1:end),'_GROUP_INDIVIDUAL.xls');
% if exist(xlsfilename, 'file') == 0   %if the excel does not exit
%     xlswrite(xlsfilename, GroupOut, 'Sheet1', 'A1');
%     warning off MATLAB:xlswrite:AddSheet
%     xlswrite(xlsfilename, GroupOutQ, 'Sheet3', 'A1');
%     xlswrite(xlsfilename, temp, 'Sheet2', 'A1');
%     xlswrite(xlsfilename, tempQ, 'Sheet4', 'A1');
%     %change sheet names
%     
%     e = actxserver ('Excel.Application'); %# open Activex server
%     %e.Visible = 1;
%     ewb = e.Workbooks.Open(xlsfilename); %# open file (enter full path!)
%     ewb.Worksheets.Item(1).Name = 'Tracking Group Mean'; %# rename sheet
%     ewb.Worksheets.Item(2).Name = 'Tracking Individual Mean'; 
%     ewb.Worksheets.Item(3).Name = 'Quant Group Mean'; 
%     ewb.Worksheets.Item(4).Name = 'Quant Individual Mean'; 
%     ewb.Save %# save to the same file
%     ewb.Close(false)
%     e.Quit
% else                 %if the excel does exit in the same folder
%     e = actxserver ('Excel.Application'); %# open Activex server
%     %e.Visible = 1;
%     ewb = e.Workbooks.Open(xlsfilename); %# open file (enter full path!)
%     ewb.Worksheets.Item(1).Cells.Clear;
%     ewb.Worksheets.Item(2).Cells.Clear;
%     ewb.Worksheets.Item(3).Cells.Clear;
%     ewb.Worksheets.Item(4).Cells.Clear;
%     ewb.Save %# save to the same file
%     ewb.Close(false)
%     e.Quit
%     xlswrite(xlsfilename, GroupOut, 'Tracking Group Mean', 'A1');
%     warning off MATLAB:xlswrite:AddSheet
%     xlswrite(xlsfilename, GroupOutQ, 'Quant Group Mean', 'A1');
%     xlswrite(xlsfilename, temp, 'Tracking Individual Mean', 'A1');
%     xlswrite(xlsfilename, tempQ, 'Quant Individual Mean', 'A1');
%     
% end
% e = actxserver ('Excel.Application'); %# open Activex server
% %e.Visible = 1;
% ewb = e.Workbooks.Open(xlsfilename); %# open file (enter full path!)
% %SheetRange1 = ewb.Sheets.Item(1).Cells;
% SheetRange2 = ewb.Sheets.Item(2).Cells;
% %SheetRange3 = ewb.Sheets.Item(3).Cells;
% SheetRange4 = ewb.Sheets.Item(4).Cells;
% SheetRange1 = ewb.Sheets.Item(1).Range('B1:L1');
% SheetRange3 = ewb.Sheets.Item(3).Range('B1:L1');
% SheetRange1.EntireColumn.HorizontalAlignment = -4108;
% SheetRange2.EntireColumn.HorizontalAlignment = -4108;
% SheetRange3.EntireColumn.HorizontalAlignment = -4108;
% SheetRange4.EntireColumn.HorizontalAlignment = -4108;
% SheetRange1.EntireColumn.AutoFit;
% SheetRange2.EntireColumn.AutoFit;
% SheetRange3.EntireColumn.AutoFit;
% SheetRange4.EntireColumn.AutoFit;
% 
% ewb.Save %# save to the same file
% ewb.Close(false)
% e.Quit


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graphs and thier associated functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% The main movement overview graph in original version, this is now an
%%% option called "plot individual VM/time colormap"
if (plotColormap)
overview_graph(fishVelocities,numWells,numFrames,usedWells,frameRate)
end

%%% The error rate graph
if (plotNoise)
    %errorRate_graph(errorPlot,usedWells,maxNoiseThresh)
    errorPlot = NOEP+TMOEP;
    allemptyWells = intersect(allfish,list(NOEP > emptyWellThresh));
    nonemptyWells = setdiff(allfish,allemptyWells);
    errorRate_graph(errorPlot,NOEP,TMOEP,usedWells,maxNoiseThresh,emptyWellThresh,nonemptyWells,allfish)
end

%%% The group mean velocity over time graph
if (plotGMVOT)
    %GMVOTplot is a [set number*frame number] matrix and each row contains 
    %the averaged frame distance for a certain fish set.
    %Meanwhile, GMVOT is a [1*frame number]matrix which contains the 
    %averaged frame distance for the last fish set
    meansmooth_graph(GMVOTplot,GMVOTsmoothplot,fishSet, 'Group Mean Instantaneous Velocity (mm/s)',frameRate);
    %figureSize = get(gcf,'Position');
    %uicontrol('Style','text',...
     %     'String','mm/s',...
     %     'Position',[10 figureSize(4)/2 50 50],...
    %      'BackgroundColor',get(gcf,'Color'),...
    %      'FontWeight','bold','FontSize',10);
    %plot(GMVOT)
end

%%% The coefficient of variation of tracking over time graph
if (plotGPTOT)
    meansmooth_graph(GPTOTplot,GPTOTsmoothplot,fishSet, 'Group Instantaneous Movement Probability (tracking)',frameRate);
end

%%% The group mean quantification over time graph
if (plotGMQOT)
    %GMQOTplot is a [set number*frame number] matrix and each row contains 
    %the averaged frame quantification data for a certain fish set.
    %Meanwhile, GMQOT is a [1*frame number]matrix which contains the 
    %averaged frame quantification data for the last fish set
    meansmooth_graph(GMQOTplot,GMQOTsmoothplot,fishSet, 'Group Mean Instantaneous Quantification (/s)',frameRate);
    %plot(GMQOT)
end

%%% The coefficient of variation of quantification over time graph
if (plotGQPTOT)
    meansmooth_graph(GQPTOTplot,GQPTOTsmoothplot,fishSet, 'Group Instantaneous Movement Probability (quantification)',frameRate);
end

%%% The intensity graphs for Vm Va and T%
if (exist('fishAreas') && plotIntensities)
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,3)), usedWells, 'Mean Velocity (mm/s)');
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,5)), usedWells, 'Active Velocity (mm/s)');
    intensity_graph(unscaledRadius, fishAreas, cell2mat(IndividOut(2:end,6)), usedWells, 'Percent Time Moving');
end

if (plotHistograms)
    histogram_graph(fishVelocities, fishSet, [1.5:.5:40],frameRate)
end

if (exist('fishCoords') && plotPPO)
    ppo_graph(fishCoords,usedWells,ppoTime,ppoTimestr)
end

%%%% Hack to generate a heatmap from the coordinates for VPConv (slow)
if (~exist('heatMap') && exist('fishCoords'))
    heatMap = zeros([max(max(fishCoords(:,:,1)))+1,max(max(fishCoords(:,:,2)))+1]);
    for i = 1:size(fishCoords,2)
        for j = 1:size(fishCoords,1)
           heatMap(floor(fishCoords(j,i,1))+1,floor(fishCoords(j,i,2))+1) = heatMap(floor(fishCoords(j,i,1)+1),floor(fishCoords(j,i,2))+1)+1;
        end
    end
end

if (exist('heatMap') && plotHeatmap)
    heatmap_graph(heatMap)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to plot overview graph (Fig 1F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function overview_graph(fishVelocities,numWells,numFrames,usedWells,frameRate)
    figure;
    imagesc(fishVelocities);
    title('Individual Velocity/Time Color Map','FontWeight','bold');
    xlabel('Time(s)');
    xlim([0 numFrames]);
    xmaxsec = numFrames/frameRate;
    set(gca,'XTick',0:300*frameRate:numFrames); 
    set(gca,'XTickLabel',0:300:xmaxsec);
    set(gca,'FontSize',8);
    ylabel('Sample');
    %set(gca,'XTickLabel', '');
    set(gca,'YTickLabel',[1,4:4:(numWells*4)]);
    set(gca,'YTick',1:4:numWells);
    hold on;
    for x = 1:numWells-1
        line([1 numFrames], [x+.5 x+.5], 'color', 'w');
    end
    colormap(flipud(gray(16)));
    plot(1,usedWells,'gs','MarkerFaceColor','g','MarkerSize',5)
    hold off;
    colorbar('YTickLabel',[0:max(max(fishVelocities))]);
    %%% User now required to manually save this image if it is wanted
    %%% saveas(gca, strcat(PathName,FileName(1:end-4),'.jpg'),'jpg');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Graph error rates (Sup 2B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function errorRate_graph(errorPlot,NOEP,TMOEP,usedWells,maxNoiseThresh,emptyWellThresh,nonemptyWells,fish)
    %%% Graph customization
    meanLineColor = [0 0 1]; %blue lines for the means and SD
    meanLineStyle = '--';
    stdLineColor = [0 0 1];
    stdLineStyle = ':';
    unusedWells = setdiff(fish,usedWells);
    
    %%%Plot the three error reports separately
    figure;  
    %total error plot
    subplot(3,3,[1 3])    
    for ii = 1:length(usedWells)
        wellcount = usedWells(ii);
        if errorPlot(wellcount) >0
            rectangle('Position',[wellcount-0.4,0,0.8,errorPlot(wellcount)],'Curvature',[0,0],...
              'FaceColor','g')
        end
        hold on
    end
    for ii = 1:length(unusedWells)
        wellcount = unusedWells(ii);
        if errorPlot(wellcount) > 0
            rectangle('Position',[wellcount-0.4,0,0.8,errorPlot(wellcount)],'Curvature',[0,0],...
              'FaceColor','r')
        end
        hold on
    end
    ylabel('% frames','FontSize',8);
    text(1,1.3*maxNoiseThresh,'Well rejection threshold','Color', 'r', 'FontSize',9);
    title('Total Tracking Errors','FontWeight','bold');
    xlim([0.5 max(fish)+0.5]);
    set(gca,'XTick',fish);
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',6);
    ylim([0 1.2*maxNoiseThresh]);   
    hold on;
    %Add the useable well threshold line
    line([0 max(fish)+1], [maxNoiseThresh maxNoiseThresh], 'color', 'r', 'LineStyle', '-'); 
    hold off;
    
    %plot the empty error report
    subplot(3,3,[4 6])
    for ii = 1:length(usedWells)
        wellcount = usedWells(ii);
        if NOEP(wellcount) > 0
            rectangle('Position',[wellcount-0.4,0,0.8,NOEP(wellcount)],'Curvature',[0,0],...
              'FaceColor','g')
        end
        hold on
    end
    for ii = 1:length(unusedWells)
        wellcount = unusedWells(ii);
        if NOEP(wellcount) > 0
            rectangle('Position',[wellcount-0.4,0,0.8,NOEP(wellcount)],'Curvature',[0,0],...
              'FaceColor','r')
        end
        hold on
    end
    ylabel('% frames','FontSize',8);
    text(1,1.3*emptyWellThresh,'Empty well threshold','Color','b','FontSize',9);
    title('No Object Errors','FontWeight','bold');
    xlim([0.5 max(fish)+0.5]);
    set(gca,'TickLength',[0 0])
    ylim([0 1.2*maxNoiseThresh]);
    set(gca,'XTick',fish);
    set(gca,'FontSize',6);    
    %Add the empty well threshold line
    hold on;
    line([0 max(fish)+1], [emptyWellThresh emptyWellThresh], 'color', 'b', 'LineStyle', '-'); 
    hold off
    
    %plot the too many objects error report
    subplot(3,3,[7 9])
    for ii = 1:length(usedWells)
        wellcount = usedWells(ii);
        if TMOEP(wellcount) > 0
            rectangle('Position',[wellcount-0.4,0,0.8,TMOEP(wellcount)],'Curvature',[0,0],...
              'FaceColor','g')
        end
        hold on
    end
    for ii = 1:length(unusedWells)
        wellcount = unusedWells(ii);
        if TMOEP(wellcount) > 0
            rectangle('Position',[wellcount-0.4,0,0.8,TMOEP(wellcount)],'Curvature',[0,0],...
              'FaceColor','r')
        end
        hold on
    end
    xlabel('Sample Number','FontSize',8);
    title('Too Many Objects Errors','FontWeight','bold');
    ylabel('% frames','FontSize',8);
    xlim([0.5 max(fish)+0.5]);
    ylim([0 1.2*maxNoiseThresh]);
    set(gca,'XTick',fish);
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',6);    
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to generate mean velocty and T% plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function meansmooth_graph(plotData, smoothplotData, fishSet, ttitle,frameRate)
    figure;
    setNum = length(fishSet);% to get the set number
    ymax = prctile(max(plotData),95);% to allow the adjustment of ylim below
    %if ymax<1 ymax=1; % stops T% graph from being at top of axis, might replace with a function to do mean +3SD or similar in future?
    %end 
    xmax = length(plotData);
    xmaxsec = xmax/frameRate;
    plotcolor = ['b','r','g','m','c','y'];
    for setCount = 1:setNum
        if setCount>6 colorCount=rem(setCount, 6); 
        else colorCount=setCount;
        end
        subplot(setNum,3,[(setCount-1)*3+1 (setCount-1)*3+3]) %adjust the subplot format automatically
        plot(plotData(setCount,:),'color',[0.8, 0.8, 0.8]);
        setCountStr = num2str(setCount);
        yaxislabel = ['Group ' setCountStr];
        ylabel(yaxislabel);
        ylim([0 1.2*ymax]);
        %text(200,ymax,yaxislabel)
        if setCount ==1 title(ttitle,'FontWeight','bold') 
        end
        if setCount == setNum
            xlabel ('Time (s)')
        end
        hold on
        plot (smoothplotData(setCount,:),'color',[plotcolor(colorCount)]);        
        xlim([0 xmax]);
        set(gca,'XTick',0:300*frameRate:xmax);
        set(gca,'XTickLabel',0:300:xmaxsec);
        set(gca,'FontSize',8);
    end
    %legend(gca,usedGroups);
    %xlim([0.5 size(usedGroups,1)+0.5]);
        %ylabel('Time (min)');
    %ylim([1 length(data)]);
    %set(gca,'YTick',floor(length(data)/10):floor(length(data)/10):length(data));
    %set(gca,'YTickLabel', floor(numFrames/2/60/10):floor(numFrames/2/60/10):floor(numFrames/2/60));
    
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function to create intensity plots (Fig 4F) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intensity_graph(unscaledRadius, fishAreas, values, usedWells, graphTitle)
    figure;
    hold on;
    %Needed for conversion from matrix to image coords. 
    maxx = max(mean(fishAreas(:,3:4)'));
    %Compute the well centers using the fish area bounding box
    fishCenters = [maxx-mean(fishAreas(:,1:2)'); mean(fishAreas(:,3:4)');];
    fishCenters = fishCenters';
    plot(fishCenters(:,2),fishCenters(:,1),'k.');

    t = linspace(0,2*pi,1000);
    r = unscaledRadius;
    j = 1;
    maxIntensity = max(values);
    minIntensity = min(values);
    for i = usedWells
        h = fishCenters(i,2);
        k = fishCenters(i,1);
        x = r*cos(t)+h;
        y = r*sin(t)+k;
        %intensity = (values(j)-minIntensity)/(maxIntensity-minIntensity);
        %intensity = values(j)/maxIntensity;
        intensity = values(j);
        set(gca,'Clim',[minIntensity,maxIntensity]);
        j = j+1;
        %fill(x,y,[intensity intensity intensity]);
        fill(x,y,intensity);
    end
    hold off;
    axis off equal;
    title(graphTitle,'FontWeight','bold');
    colormap(flipud(gray(128)));
    %colorbar('YTickLabel',floor([minIntensity:floor(maxIntensity-minIntensity)/8:maxIntensity].*100)/100);
    colorbar('YTickLabel',[minIntensity:((maxIntensity-minIntensity)/8):maxIntensity]);
end

function histogram_graph(fishVelocities, fishSet, range, frameRate)    
    setNum = length(fishSet);
    ylimtem = zeros(1,setNum);
    figure
    for setCount = 1:setNum
        startdot = 1;
        AxesHandle(setCount) = subplot(setNum,3,[(setCount-1)*3+1 (setCount-1)*3+3]);
        setWell = fishSet{setCount};
        instV = reshape(fishVelocities(setWell,:), 1, size(fishVelocities,2)*length(setWell)); 
        instV = sort(instV);
        histDist = histc(instV,range);    
        while instV(startdot) < range(1);
             startdot = startdot+1;
        end
        startdot = startdot-1;
        count = zeros(1,length(range)+1);
        wHistDist = zeros(1,length(range));
        count(1) = startdot;
        for i = 1:length(range)
             count(i+1) = count(i)+histDist(i);
             wHistDist(i) = sum(instV(count(i):count(i+1)))/frameRate; %calculate the distance for each velocity range, in (mm)
        end       
        wHistDist = wHistDist./length(setWell); %notmalize to one fish
        bar(range,wHistDist);
        xlim([min(range) max(range)]);
        if setCount ==1 title('Total Displacement (mm)/Instantaneous Velocity (mm/s)','FontWeight','bold') 
        end
        if setCount == setNum
            xlabel ('Velocity (mm/s)')
        end
        setCountStr = num2str(setCount);
        yaxislabel = ['Group ' setCountStr];
        ylabel(yaxislabel)
        allYLim = get(AxesHandle, {'YLim'});
        allYLim = cat(2, allYLim{:});
    end
    set(AxesHandle, 'YLim', [min(allYLim), 1.2*max(allYLim)]);
    %ylim([0 1.2*max(ylimtem)])
    %%% Allow matlab to autoscale
    %ylim([0 3000]);
end

function ppo_graph(fishCoords, wells, range, rangestr)
    figure;  
    title(strcat('Vectors for Frames -',rangestr),'FontWeight','bold');
    if (range(1)>0 && range(end) < size(fishCoords,2))
        plotPathOverlay(fishCoords(wells,range,:));
    elseif (range(1) && range(end) >= size(fishCoords,2))
        range = range(1):size(fishCoords,2);
        plotPathOverlay(fishCoords(wells,range,:));
    else
        plotPathOverlay(fishCoords(wells,:,:));
    end
end

function heatmap_graph(heatMap)
    figure;
    %colormap(hot(128));
    lineheatMap = reshape(heatMap, 1, size(heatMap,1)*size(heatMap,2)); 
    percent5 = prctile(lineheatMap,5);
    percent95 = prctile(lineheatMap,95);
    colormap(flipud(gray));
    image(heatMap);
    caxis ([percent5 percent95])
    colorbar
    daspect([1 1 1]);
    axis off;
    title('Positional Heatmap','FontWeight','bold');
end

function standerror = ste(individData)
    standerror = nanstd(individData)/sqrt(length(individData));
end


    
