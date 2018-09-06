function [m,n, stepX, stepY, possibleConfs] = matrixLayoutPointsInABox(nbPoints, sizeBoxX, sizeBoxY)
%How to best split a number of points in the area of a square,
%boundaries being excluded, IN A MATRIX ARRANGED WAY.
%askiitians.com/forums/Algebra/22/34781/matrices.htm

%For a RANDOM ARRANGED WAY :
%stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly/50746409#50746409

% %To better understand with a non-empty workspace:
% %example : comment functiontitle line & the 'end function'
% %uncomment the following line, press f5.
% nbPoints = 46 ; sizeBoxX = 1 ; sizeBoxY =0.9 ;

divisors = [];
for i=1:nbPoints
    if mod(nbPoints,i) == 0 ; divisors(end+1) = i ; end
end % hence 

possibleConfs = []; % all [m,n] possible.
if mod( numel(divisors) , 2 ) == 0 %length is even (pair in french)
    while numel(divisors) ~= 0
        possibleConfs(end+1,:) = [divisors(1) divisors(end)];
        possibleConfs(end+1,:) = [divisors(end) divisors(1)];
        divisors(1) = [] ; divisors(end) = [];
    end
else %length is odd (impair in french)
    while numel(divisors) ~= 0
        if numel(divisors) == 1
            possibleConfs(end+1,:) = [divisors(1) divisors(1)];
            divisors(1) = [];
        else
            possibleConfs(end+1,:) = [divisors(1) divisors(end)];
            possibleConfs(end+1,:) = [divisors(end) divisors(1)];
            divisors(1) = [] ; divisors(end) = [];
        end
    end
end
    
steps = [sizeBoxX./(possibleConfs(:,2)+1) sizeBoxY./(possibleConfs(:,1)+1)] ; %[dx dy]

%I take the min dx or dy for each doublet, and select the doublet
%containing the max of these mins :
[~, selectedRow] = max( min(steps') );

selectedConf = possibleConfs(selectedRow,:);
m = selectedConf(1) ; n = selectedConf(2);

selectedSteps = steps(selectedRow,:);
stepX = selectedSteps(1) ; stepY = selectedSteps(2);

end