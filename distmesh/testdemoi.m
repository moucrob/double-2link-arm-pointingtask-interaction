clc ; %clf(2)
% digitsOld = digits(100);
nbPoints = 46;
% factor = 100000;
boxX = 1; %boxX = boxX*factor;
boxY = 1; %boxY = boxY*factor;
boundingbox = [0,0;boxX,boxY]; %xmin ymin ; xmax ymax
verticesWeWant = [0,0;boxX,0;boxX,boxY;0,boxY];
fd=@(p) drectangle(p,0,boxX,0,boxY); %xmin xmax ; ymin ymax

flag1 = 0; flag2 = 0; %temporary, used if script is a script and not a function
step = 0.17666666666668;
decrease = 0.01;
% step = 20000;
decrease% = 100;
decimals = 100;

tic;

figure(1)
[p,t]=distmesh2d(fd,@huniform,step,boundingbox,verticesWeWant);
disp(['size(p,1) = ',num2str(size(p,1))])
if size(p,1) < nbPoints
    
    sizes = [inf]; stepAboveNbPoints = 0; %just to enter the loop
    stepUnderNbPoints = step; xMin = 0; xMax = nbPoints+100;
    while any(unique(sizes) ~= nbPoints)
        sizes = [];
        while 1
            step = step - decrease;
            figure(1)
            [p,t]=distmesh2d(fd,@huniform,step,boundingbox,verticesWeWant);
            disp(['size(p,1) = ',num2str(size(p,1))])
            if size(p,1) == nbPoints
                flag1 = 1 %we have placed our desired nbPoints
                pSave = p;
            end
            disp(['stepAboveNbPoints = ',num2str(stepAboveNbPoints)])
            disp(['stepUnderNbPoints = ',num2str(stepUnderNbPoints)])
            figure(2) ; grid minor ; set(gca, 'YScale', 'log')
            plot(size(p,1),step,'.'); hold on;
            Limits = [stepAboveNbPoints stepUnderNbPoints];
            ylim(Limits)
            xlim([xMin xMax]) ; drawnow
            sizes(end+1) = size(p,1);
            if size(p,1) > nbPoints
                xMax = size(p,1);
                disp('cycle')
                %use directly return if script is a function, it will
                %serves as double break !
                stepAboveNbPoints = step;
                flag2 = 1; %we have found the most optimal way of placing our nbPoints
                pAtLeastIHaveTried = p;
                break
            end
            xMin = size(p,1);
        end
        if (flag1==1) && (flag2 == 1) ; break ; end
        stepUnderNbPoints = stepAboveNbPoints + decrease;
        decrease = (stepUnderNbPoints-stepAboveNbPoints)/2;
        step = stepUnderNbPoints;
    end
    
end

tEnd = toc;
disp(['tEnd = ',num2str(tEnd)])