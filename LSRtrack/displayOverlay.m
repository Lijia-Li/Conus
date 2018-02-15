function displayOverlay(Img,bwImage,target)
    % Create a 3 channel output image based on the full grayscale input
    % frame
    %bwImage = im2bw(bwImage);
    %Img = Img./2;
    redOut = Img;
    greenOut = Img;
    blueOut = Img;
    % To the green channel, increase intensity for all seach space
    % greenOut(not(bwImage)) = greenOut(not(bwImage)) + max(max(Img))/10;
    % To the blue channel, increase intensity for all target objects
    blueOut(target) = blueOut(target)./2 + max(max(blueOut(target)));
    % To the red channel, increase intensity for thresholded out pixels
    redOut(bwImage) = redOut(bwImage)./2 + max(max(redOut(bwImage)));
    % Create the image with overlay
    alphaImage = cat(3,redOut,greenOut,blueOut);
    
    alphaImage = uint8(alphaImage);
    if (max(max(max(alphaImage)))<255)
        alphaImage = alphaImage.*127;
    end
    % Display it, and draw cenroid marks in blue
    %maxValue = max(max(max(alphaImage)));
    %alphaImage = alphaImage ./ maxValue;
    %keyboard;
    imagesc(uint8(alphaImage));
end
