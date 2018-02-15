function runTracking3(handles)
%%%% VERSION 3.99 6/21/11
%%%% Optimized Windows/Mac (Also works with some Unix systems)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get input parameters from GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%watchFlag = get(handles.displayChk,'Value');
%watchWell = str2double(get(handles.watchWell,'String'));
% scaleFactor: The percent of the total area that is used for tracking
scaleFactor = str2double(get(handles.scaleFactor,'String'));
% wellThresh: The number of grouped pixels needed to be considered a well 
wellThresh = str2double(get(handles.wellThresh,'String'));
% fishThresh: The number of grouped pixels needed to be a fish
fishThresh = str2double(get(handles.fishThresh,'String'));
% minimumMovement: The smallest recordable fish movement
minimumMovement = str2double(get(handles.minimumMovement,'String'));
% trackingThresh: The pixel grey level cutoff for a fish
trackingThresh = str2double(get(handles.trackingThresh,'String'));
% quantingThresh: The grey scale cut off for a fish in quantization
quantingThresh = str2double(get(handles.quantingThresh,'String'));
% fileList: The list of movies to track
fileList = get(handles.fileList,'String');
% directoryName: The directory the tracking videos are stored
directoryName =  get(handles.directory,'String');
% outputPath: Where output information will be written
outputPath = get(handles.OutputPath,'String');
% alignmentFreq: How often well coordinates are updated
alignmentFreq = str2double(get(handles.alignFreq,'String'))*.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check Input and Output Information 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numMovies = 0;
if (isempty(fileList) || isempty(directoryName) || strcmp(outputPath,'No directory selected...') || strcmp(get(handles.FilesPath,'String'),'No files selected...'))
    errordlg('You did not specify input file(s) or an output directory');
    set(handles.goBut,'Enable','on');
    return;
end
if (iscell(fileList))
    files = cat(1, char(fileList(:)));
    numMovies = length(fileList);
else
    files = fileList;
    numMovies = 1;
end

set(handles.status,'String','Movie information is being read. Please wait....');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Main Loop (for each movie)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for currentMovie = 1:numMovies
    
    %----- Step 1: Load frame 1 and threshold to black and white
    set(handles.CurrentMovie,'String',strcat(int2str(currentMovie),'/',int2str(numMovies)));
    %Open movie using the mplayer function (not Matlab built in)
	 %{
    [av_hdl, av_inf] = mplayerOpen([directoryName, files(currentMovie,:)]);
    if ~isempty(av_hdl)
            %Get the first frame
			firstFrame = mplayerReadMex(av_hdl, 1);
            numFrames = av_inf.NumFrames;
            frameRate = av_inf.fps;
            currentFrame = reshape(firstFrame/255,[av_inf.Height,av_inf.Width,3]);
    else
            fprintf('Could not open movie!');
            return;
    end
	 %}
    readerobj = VideoReader(strcat(directoryName, files(currentMovie,:)));
    numFrames = readerobj.NumberOfFrames;
    frameRate = readerobj.FrameRate;
    currentFrame = read(readerobj,1);
    currentFrame = rgb2gray(currentFrame);
    lastFrame = currentFrame;
    bwFrame = im2bw(currentFrame);
    bwFrame = bwareaopen(bwFrame,wellThresh);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Align the background (update well coordinates)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %----- Step 2-8: Get well coordinates
    [unscaledRadius,radius,fishAreas,background] = alignBackground(handles,bwFrame, scaleFactor);
    if (fishAreas == -1)
        return;
    end
    numWells = size(fishAreas,1);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Tracking Loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare output data structures 
    % fishCoords contains the xy coords for each fish per frame
    fishCoords(numWells,numFrames,2) = 0;
    % fishDistances contains the movement since the last frame for each
    % fish
    fishDistances(numWells,numFrames,1) = 0;
    % fishTotalDistance contains the cummulative distance for each fish
    fishTotalDistance(numWells) = 0;
    % errorCount contains the number of encountered bad frames
    errorCount = 0;
    % relativeFishLoc contains the fishes position relative to the well
    relativeFishLoc(numWells,:) = [0,0];
    % noObjectError contains the number of times the fist was lost in
    % tracking for each well
    noObjectError(numWells) = 0;
    noObjectErrorByFrame(numWells,numFrames)=0; % Line 196 increments for each fish
    
    % tooManyObjectError contains the number of times more than one object
    % was found in each well
    tooManyObjectError(numWells) = 0;
    tooManyObjectErrorByFrame(numWells,numFrames)=0;
    % heatMap, an image of where the fish spend their time
    heatMap(size(bwFrame,1),size(bwFrame,2)) = 0;
    totalQuant(size(bwFrame,1),size(bwFrame,2)) = 0;

    skip = 1; % Don't skip any frames
    for frameNum = skip:skip:(numFrames-1)
        oldFN = frameNum;
        frameNum = frameNum/skip;
        frameTime = tic();
        %Allow for real time variable changes
        fishThresh = str2double(get(handles.fishThresh,'String'));
        minimumMovement = str2double(get(handles.minimumMovement,'String'));
        trackingThresh = str2double(get(handles.trackingThresh,'String'));
        quantingThresh = str2double(get(handles.quantingThresh,'String'));
        alignmentFreq = str2double(get(handles.alignFreq,'String'))*.01;
            
        set(handles.status,'String','Tracking...');
        set(handles.CurrentFrame,'String',strcat(int2str(frameNum),'/',int2str(numFrames)));

        %----- Step 9: Load the next frame, threshold to black and white 
        try 
            %currentFrame = mplayerReadMex(av_hdl, frameNum);
				currentFrame = read(readerobj,frameNum);
        catch ME
            set(handles.status,'String',sprintf('%s %d','Error reading frame: ',frameNum));
            errorCount = errorCount+1;
            if (errorCount > 5)
                warndlg('Too many bad frames in this video');
                return;
            else
                continue;
            end
        end
        %currentFrame = reshape(currentFrame/255,[av_inf.Height,av_inf.Width,3]);
        grayFrame = rgb2gray(currentFrame);
        %count the pixels which become brighter
        quant = grayFrame-lastFrame;
        quant(quant<quantingThresh) = 0;
        %totalQuant = totalQuant+quant;
        lastFrame = grayFrame;
        %grayFrame = grayFrame .* 255;
        bwFrame = im2bw(grayFrame,trackingThresh);
        bwFrame = bwareaopen(bwFrame,wellThresh);
        
        %---- Step 15: Plate Alignment rescheduled? (moved to use newest bwFrame)
        % Align the plate again if so
        if ((mod(frameNum,uint32(alignmentFreq*numFrames)) == 0) || strcmp('True',get(handles.reAlign,'String')))
            
            scaleFactor = str2double(get(handles.scaleFactor,'String'));
            wellThresh = str2double(get(handles.wellThresh,'String'));
            set(handles.reAlign,'String','False');
            
            %Using automatic threshold is better for finding background....
            [unscaledRadius,radius,fishAreas,background] = alignBackground(handles,bwFrame, scaleFactor);
            
            if (fishAreas == -1)
                return;
            end
            numWells = size(fishAreas,1);
            save(strcat(outputPath,'/',files(currentMovie,:),'.mat'),'fishDistances', 'fishCoords', 'noObjectError', 'fishAreas','frameRate','radius', 'tooManyObjectError','heatMap','unscaledRadius', 'noObjectErrorByFrame','tooManyObjectErrorByFrame');
        end
        
        bwFrame = not(bwFrame);
        %----- Step 10: Subtract plate image for both tracking and quant
        bwFrame = bwFrame.*background;
        quant = double(quant).*background;
        %----- Step 11: Remove small artifacts
        bwFrame = bwareaopen(bwFrame,fishThresh);
        %%%%%%%
        

        %----- Step 12: Locate fish in each well
        for wellNum = 1:numWells      
            % currentWellLoccontains the current fishAreas coords (Row1,Row2,Col1,Col2)
            currentWellLoc = fishAreas(wellNum,:);
            %fish contains the fish pixels in the fishAreas
            fish = bwFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
            fishQ = quant(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
            %find the largest object (fish) by area, store the centroid (relative)
            % if no objects are found, continue
            wellObject = regionprops(fish,'Centroid','Area');
            [unused, order] = sort([wellObject(:).Area],'descend');
            wellObject = wellObject(order);

            %----- Step 13: Does the number of 'larvae' in each well = 1?
            %----- Larval count = 0
            if (isempty(wellObject))
		        noObjectErrorByFrame(wellNum,frameNum)=1;
                noObjectError(wellNum) = noObjectError(wellNum) + 1;
                if (frameNum == 1) %if first frame, use well center, otherwise use the last fish coords for this well
                    wellObject(1).Centroid(2) = round((currentWellLoc(2)-currentWellLoc(1))/2);
                    wellObject(1).Centroid(1) = round((currentWellLoc(4)-currentWellLoc(3))/2);
                else
                    wellObject(1).Centroid(2) = relativeFishLoc(wellNum,1);
                    wellObject(1).Centroid(1) = relativeFishLoc(wellNum,2);
                end
                %{
                %%%%%%% ENABLE FOR MANUAL FISH SELECTION
                fish =[27,28,29,30,35,36,37,38,43,44,45,46,51,52,53,54,59,60,61,62,67,68,69,70;];
                if ismember(wellNum,fish)
                    set(get(handles.WatchWellFig,'parent'),'CurrentAxes',handles.WatchWellFig);
                    set(handles.watchWell,'String',int2str(wellNum));
                    currentWellLoc= fishAreas(wellNum,:);
                    target = bwFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
                    displayOverlay(grayFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4)),not(background(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4))), target);
                    axis off square;
                    hold on;    
                    plot(relativeFishLoc(wellNum,2), relativeFishLoc(wellNum,1), 'g+');
                    text(2,4, strcat('No objects detected in well ',num2str(wellNum)), 'BackgroundColor', [.7 .9 .7]);
                    hold off;
                    drawnow;
                    pos = ginput(1);
                    deltaX = abs(relativeFishLoc(wellNum,2)) - pos(1);
                    deltaY = abs(relativeFishLoc(wellNum,1)) - pos(2);
                    deltaDist = (sqrt(deltaX^2 + deltaY^2));
                    wellObject(1).Centroid(1) = pos(1);
                    wellObject(1).Centroid(2) = pos(2);

                    %%%% Well, Frame, Dist, X, Y
                    LostObjectError(1,:) = [0,0,0,0,0];
                    LostObjectError(end+1,:) = [wellNum,frameNum,deltaDist,pos(1),pos(2)];
                end
                %%%%%%%
                %}
            %----- Larval count > 1
            elseif (length(wellObject) > 1)
                tooManyObjectErrorByFrame(wellNum,frameNum)=1;
                tooManyObjectError(wellNum) = tooManyObjectError(wellNum) + 1;
                %{
                %%%%%%% Code to verify the right object was selected as the
                %%%%%%% fish
                fish = [27,28,29,30,35,36,37,38,43,44,45,46,51,52,53,54,59,60,61,62,67,68,69,70;];
                if ismember(wellNum,fish)           
                    relativeFishLoc(wellNum,:) = [wellObject(1).Centroid(2), wellObject(1).Centroid(1)];
                    set(get(handles.WatchWellFig,'parent'),'CurrentAxes',handles.WatchWellFig);
                    set(handles.watchWell,'String',int2str(wellNum));
                    currentWellLoc = fishAreas(wellNum,:);
                    target = bwFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
                    displayOverlay(grayFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4)),not(background(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4))), target);
                    axis off square;
                    hold on;    
                    plot(relativeFishLoc(wellNum,2), relativeFishLoc(wellNum,1), 'g+');
                    text(2,4, strcat('Too many objects detected in well ',num2str(wellNum)), 'BackgroundColor', [.7 .9 .7]);
                    hold off;
                    drawnow;
                    pause;
                end
                %%%%%%%
                %}
            end
            % Use the largest objects 
            relativeFishLoc(wellNum,:) = [wellObject(1).Centroid(2), wellObject(1).Centroid(1)];

            %----- Step 14: Store the absolute coords of the current fish for this frame
            Col = relativeFishLoc(wellNum,2)+currentWellLoc(3);
            Row = relativeFishLoc(wellNum,1)+currentWellLoc(1);
            fishCoords(wellNum,frameNum,:) = [Row,Col];
            heatMap(floor(Row),floor(Col)) = heatMap(floor(Row),floor(Col)) + 1;
            % Compute distance traveled since last frame and store
            if (frameNum>1)
                deltaX = abs(fishCoords(wellNum,frameNum,2) - fishCoords(wellNum,frameNum-1,2));
                deltaY = abs(fishCoords(wellNum,frameNum,1) - fishCoords(wellNum,frameNum-1,1));
                deltaDist = (sqrt(deltaX^2 + deltaY^2));
                if (deltaDist > minimumMovement)
                    fishDistances(wellNum,frameNum) = deltaDist;
                    fishTotalDistance(wellNum) = fishTotalDistance(wellNum) + deltaDist;
                else
                    fishDistances(wellNum,frameNum) = 0;
                end
            end
            % Compute the quant value for this fish
            fishQuants(wellNum,frameNum) = sum(sum(fishQ));
            
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Display the tracking in one well if requested
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (get(handles.displayChk,'Value') == 1)
            figure(get(handles.WatchWellFig,'parent'));
            set(get(handles.WatchWellFig,'parent'),'CurrentAxes',handles.WatchWellFig);
            watchWell = str2double(get(handles.watchWell,'String'));
            if (watchWell > numWells)
                watchWell = 1;
                set(handles.watchWell,'String',int2str(watchWell));
            end
            if (watchWell < 1)
                watchWell = numWells;
                set(handles.watchWell,'String',int2str(watchWell)); 
            end
            currentWellLoc= fishAreas(watchWell,:);
            target = bwFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
            % displayOverlay is given the image gray(Col1=>Col2,Row1=>Row2) [grey], 
            %  the background mask from (Col1=>Col2,Row1=>Row2) [red],
            %  and the target (Col1=>Col2,Row1=>Row2) [blue]
            if (get(handles.wTrack,'Value'))
                displayOverlay(grayFrame(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4)), not(background(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4))), target);
                axis off square;
                hold on;    
                % Plot the fish's centroid with a red dot if it's below the moving threshold and a yellow dot if it isn't
                if (fishDistances(watchWell,frameNum) > minimumMovement)
                    plot(relativeFishLoc(watchWell,2), relativeFishLoc(watchWell,1), 'o','MarkerEdgeColor','k',...
                        'MarkerFaceColor','g',...
                        'MarkerSize',6);
                else
                    plot(relativeFishLoc(watchWell,2), relativeFishLoc(watchWell,1), 'o','MarkerEdgeColor','k',...
                        'MarkerFaceColor','r',...
                        'MarkerSize',6);
                end
            elseif (get(handles.wQuant,'Value'))
                axis off square;
                hold on; 
                
                %quant = 1-(quant/max(max(quant)));
                %target(:,:)=0;
                %displayOverlay(quant(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4)), target, not(background(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4))));
                
                wellquant = quant(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
                colormap(gray(128))
                %colorbar
                reversedquant = 128-wellquant;
                reversedquant = reversedquant.*background(currentWellLoc(1):currentWellLoc(2),currentWellLoc(3):currentWellLoc(4));
                imagesc(reversedquant);
                caxis([0 128])
                
            end

            % Display the fish's distances:
            text(2,size(target,2)*.05, sprintf('Last Distance: %.2f',fishDistances(watchWell,frameNum)), 'BackgroundColor', [.7 .9 .7]);
            text(2,size(target,2)*.1, sprintf('Total Distance: %.2f',fishTotalDistance(watchWell)), 'BackgroundColor', [.7 .9 .7]);
            %text(size(target,1)*.75,size(target,2)*.05, sprintf('Last Quant: %.2f',fishQuants(watchWell,frameNum)), 'BackgroundColor', [.7 .9 .7]);
            %text(size(target,1)*.75,size(target,2)*.1, sprintf('Total Quant: %.2f',sum(fishQuants(watchWell,:))), 'BackgroundColor', [.7 .9 .7]);
            text(size(target,1)*.65,size(target,2)*.05, sprintf('Last Location: [%.2f,%.2f]',fishCoords(watchWell,frameNum,2), fishCoords(watchWell,frameNum,1)), 'BackgroundColor', [.7 .9 .7]);
            
            hold off;
            drawnow;
        end

        %Compute fps and remaining time, display in window
        fps = toc(frameTime);
        fps = 1/fps;
        set(handles.fps,'String',num2str(fps));
        timeLeft = ((numFrames-frameNum)/fps)/60;
        set(handles.RemainingTime,'String',num2str(timeLeft));
        drawnow;
    end
    set(handles.CurrentFrame,'String',strcat(int2str(numFrames),'/',int2str(numFrames)));
    set(handles.RemainingTime,'String','0.0');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Write the output files
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dlmwrite(strcat(outputPath,'/',files(currentMovie,:),'.dist'), fishDistances','newline','unix');
    imwrite(background,strcat(outputPath,'/',files(currentMovie,:),'.jpg'),'jpg');
    save(strcat(outputPath,'/',files(currentMovie,:),'.mat'),'fishDistances', 'fishCoords', 'fishQuants', 'noObjectError', 'fishAreas','frameRate','radius', 'tooManyObjectError','heatMap','unscaledRadius', 'noObjectErrorByFrame','tooManyObjectErrorByFrame');
    %clear vars not needed for the next video and loop
    clear background fishAreas fishCoords fishDistances fishQuants fishTotalDistance ;
end
h = msgbox('Video(s) have finished tracking!','Tracking Completed');
uiwait(h);
end

% Input: the black and white frame, the well area scale factor
%Output: The unscaled radius, radius, fish tracking areas, and background
function [unscaledRadius,radius,fishAreas,background] = alignBackground(handles, bwFrame, scaleFactor)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Find the wells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(handles.status,'String','Aligning well locations....');
    frameHeight = size(bwFrame,1);
    frameWidth = size(bwFrame,2);
    
    %----- Step 2: Locate potential 'wells'
    Well = regionprops(bwFrame,'Centroid','Area','BoundingBox');
    numWells = length(Well);
    %----- Step 3: Find median 'well' area 
    % Wells should dominate image, so the median object area should
    % correspond to them
    medianWellArea = median([Well(:).Area]);
    %----- Step4: Select wells (within 20% of the median 'well' area
     upperWellArea = medianWellArea * 1.2;
     lowerWellArea = medianWellArea * 0.2;
     selected = [];
     for i = 1:numWells 
         if (Well(i).Area > lowerWellArea && Well(i).Area < upperWellArea)
             selected(end+1) = i;
         end
     end
    Well = Well(selected);
    numWells = length(Well);
    
    %----- Step 5: Is the number in the set 3*2^n?
    %This number corresponds to typical plate arrangements: 6,24,48,96 well
   % if ~ismember(numWells,3*2.^[1:10])
    %    warndlg(strcat('Irregular Well number. Found: ',int2str(numWells)));
     %   figure;
      %  colormap(gray(128));
       % bwFrame = 128*bwFrame;
        %image(bwFrame);
  %      axis off equal;
   %     hold on;
    %    radius = sqrt(medianWellArea/pi);
     %   for j = 1:numWells
      %      circle(Well(j).Centroid,radius,1000);
       %     text(Well(j).Centroid(1),Well(j).Centroid(2), int2str(j), 'FontSize',8, 'HorizontalAlignment', 'Center');
        %end   
  %      hold off;
   %     drawnow;
    %    fishAreas = -1;
     %   background = -1;
      %  return;
   % end
    
    %----- Step 6: Compute mean well area and radius
    meanWellArea = mean([Well(1:round(numWells/2)).Area])*scaleFactor;
    radius = sqrt((meanWellArea*scaleFactor)/pi);
    unscaledRadius = sqrt((meanWellArea)/pi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sort & Show: sort by row, then column. Allow user to verify in GUI
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:numWells
        Col(i) = Well(i).Centroid(1);
        Row(i) = Well(i).Centroid(2);
    end
    [unused, order] = sort(Col);
    Row = Row(order);
    Well = Well(order);
    % Column sort
    numWellRows = floor(sqrt((numWells*2)/3));
    % Correct for a couple plate well arragements 
    if (numWells == 12)
        numWellRows = 3;
    elseif (numWells == 48)
        numWellRows = 6;
    end
    
    for i = 1:numWellRows:numWells 
        [unused, order] = sort(Row(i:i+numWellRows-1));
        SortedWell(i:i+numWellRows-1,:) = Well(order+(i-1));
    end
    Well = SortedWell;
    
    set(get(handles.AlignAxisFig,'parent'),'CurrentAxes',handles.AlignAxisFig);
    colormap(gray(128));
    bwFrame = 128*bwFrame;
    image(bwFrame);
    hold on;
    for j = 1:size(SortedWell)
        circle(Well(j).Centroid,radius,1000);
        text(Well(j).Centroid(1),Well(j).Centroid(2), int2str(j), 'FontSize',8, 'HorizontalAlignment', 'Center');
    end   
    hold off;
    axis off equal;
    drawnow;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Search Areas: generate a search area for each well
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fishAreas contains the four corners of the bounding box [Col1,Row1,dCol+Col1,dRow+Row1] for each fish (numWell)
    fishAreas(numWells,:) = [0,0,0,0];
    % background is the image to subtract from tracking boxes
    background(frameHeight, frameWidth) = 0; 
    for wellNum = 1:numWells
        % Coords for searching
        Col1 = uint32(Well(wellNum).BoundingBox(1));
        Row1 = uint32(Well(wellNum).BoundingBox(2));
        Col2 = uint32(Well(wellNum).BoundingBox(1)+Well(wellNum).BoundingBox(3));
        Row2 = uint32(Well(wellNum).BoundingBox(2)+Well(wellNum).BoundingBox(4));
        for i = Col1:Col2
            for j = Row1:Row2
                if inCircle(Well(wellNum).Centroid,radius,[i,j])
                    % Pixel is inside the search region
                    background(j,i) = 1;
                end
            end
        end
        % Remember coords for each fish well
        fishAreas(wellNum,:) = [Row1,Row2,Col1,Col2];
    end
end

% Input: the center(x,y) and radius of a circle, and a point(x,y)
%Output: 1 if the point falls within the circle, 0 otherwise
function state = inCircle(center,radius,point)
    xDist = double(abs(double(point(1))-center(1)));
    yDist = double(abs(double(point(2))-center(2)));
    distance = sqrt(xDist^2+yDist^2);
    if (distance > radius)
        state = 0;
    else
        state = 1;
    end
end


function H=circle(center,radius,NOP,style)
%--------------------------------------------------------------------------
% H=CIRCLE(CENTER,RADIUS,NOP,STYLE)
% This routine draws a circle with center defined as
% a vector CENTER, radius as a scaler RADIS. NOP is 
% the number of points on the circle. As to STYLE,
% use it the same way as you use the rountine PLOT.
% Since the handle of the object is returned, you
% use routine SET to get the best result.
%
%   Usage Examples,
%
%   circle([1,3],3,1000,':'); 
%   circle([2,4],2,1000,'--');
%
%   Zhenhai Wang <zhenhai@ieee.org>
%   Version 1.00
%   December, 2002
%--------------------------------------------------------------------------
if (nargin <3),
 error('Please see help for INPUT DATA.');
elseif (nargin==3)
    style='b-';
end;
THETA=linspace(0,2*pi,NOP);
RHO=ones(1,NOP)*radius;
[X,Y] = pol2cart(THETA,RHO);
X=X+center(1);
Y=Y+center(2);
H=plot(X,Y,style);
axis square;
end

