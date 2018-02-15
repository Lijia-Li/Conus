function smoothedData = gaussSmooth(data,wSize)
%%%% VERSION 1.2 12/02/10  
%%%% Windows/Mac/Unix
%%%% This function is used to apply a gaussian filter to a 2D array
%%%% The input data is an arrary with more columns than rows (ie. 96x32000)
%%%% The input wSize is the integer size of the smoothing window 
%%%% (odd numbers work best). The larger the window, the more smoothing. 
%%%% NOTE: this function depends on signal processing toolbox (gausswin),
%%%% though this function can be easily written if you don't have it. 
%convert to odd size if neccessary
	if (mod(wSize,2)==0)
	    wSize = wSize+1;
    end
	%%% Initialize variables
	halfSize = (wSize-1)/2;
	gaussCoeffs = gausswin(wSize);
	numRows = size(data,1);
	numCols = size(data,2);
	smoothedData(numRows,numCols) = 0;

	for pos = 1:numCols
	    if (pos<=halfSize)
		offset = halfSize-pos+1;
		scalr = 1/sum(gaussCoeffs(1+offset:wSize));
		filter = gaussCoeffs(1+offset:wSize)*scalr;
		smoothedData(:,pos) = data(:,1:wSize-offset)*filter;
	    elseif (pos>numCols-halfSize)
		offset = pos-(numCols-halfSize);
		scalr = 1/sum(gaussCoeffs(1:wSize-offset));
		filter = gaussCoeffs(1:wSize-offset)*scalr;
		smoothedData(:,pos) = data(:,pos-halfSize:end)*filter;
	    else
		scalr = 1/sum(gaussCoeffs);
		smoothedData(:,pos) = data(:,pos-halfSize:pos+halfSize)*gaussCoeffs.*scalr;
	    end
        end

end
