function angVelocities = plotPathOverlay(coords, backgroundPath)
%%%% VERSION 1.2 12/02/10  
%%%% Windows/Mac/Unix
%%%% This function is used to display the vectors given by successive x,y 
%%%% coordinates in coords, and optionally displays them over an image
%%%% specified by the location backgroundPath. This function is time
%%%% intensive. coords should be an array like: [sample#,frame#, [x,y]] 
%%%% backgroundPath can optionally given. 
	if (nargin == 2)
	    background = imread(backgroundPath,'jpg');
	    imagesc(background);
	    colormap(cool(2));
	    hold on;
	end
	%Size returns [#samples, #frames, 2]
	angVelocities(size(coords,1),size(coords,2)) = 0;
	%numRows = sqrt(size(coords,1)*2/3);
	%numCols = size(coords,1)/numRows;
	%minX = min(min(coords(:,:,1)));
	%maxX = max(max(coords(:,:,1)));
	%minY = min(min(coords(:,:,2)));
	%maxY = max(max(coords(:,:,2)));
	%TextOffset = (maxX-minX)/numCols;
	%yTextOffset = (maxY-minY)/numRows;
	%xTextCoord = minX;
	%yTextCoord = minY;
	hold on;
	%figure;
	for i = 2:size(coords,2)-1
	    for sample = 1:size(coords,1)
            disp(coords)
            disp(i)
            disp(sample)
            x1 = coords(sample,i-1,2);
            x2 = coords(sample,i,2);
            y1 = coords(sample,i-1,1);
            y2 = coords(sample,i,1);
            if (x1~=0 && x2~=0 && y1~=0 && y2~=0)
                line([coords(sample,i-1,2),coords(sample,i,2)],[coords(sample,i-1,1),coords(sample,i,1)],'Color','k');	    
                Vect1 = [coords(sample,i,1) - coords(sample,i-1,1), coords(sample,i,2) - coords(sample,i-1,2), 0];
                Vect2 = [coords(sample,i+1,1) - coords(sample,i,1), coords(sample,i+1,2) - coords(sample, i,2), 0];
                %dotProd = dot(Vect1,Vect2);
                %magProd = norm(Vect1)*norm(Vect2);
                %returns NaN if ans = 0, use atan version instead
                %angVel = acos(dotProd/magProd);
                angVel = atan2(norm(cross(Vect1,Vect2)),dot(Vect1,Vect2));
                angVelocities(sample,i) = (pi*angVel)/180;
            end
	    end
	end

	%{
	for sample = 1:size(coords,1)
	    text(mean(coords(sample,:,2)),mean(coords(sample,:,1)), int2str(sample));
	end
	%}
	axis off equal;
	axis ij
end
