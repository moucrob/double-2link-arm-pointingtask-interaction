%% Credits to:
%user104676 :
%https://stackoverflow.com/questions/401847/circle-rectangle-collision-detection-intersection
%However he distincts case where 
%-circle intersect any edge?
%-any vertex penetrates the circle?
%Where I merge those two cases thanks to the function
%do1Segment1CircleIntersect
%% Workspace for verifying the function (comment the function title and ending lines):
% VerticesNx2 = [2 1 ; 1 2 ; 1 3 ; 4 3];
% r = 0.05;
% nbIntervalsForCircle = 50;
% sentence = {'No collision','COLLISION'};
% 
% set(0,'DefaultFigureWindowStyle','docked')
% figure;
% set(0,'DefaultFigureWindowStyle','normal')
% patch(VerticesNx2(:,1),VerticesNx2(:,2),'b')
% axis square ; xlim([-1 6]) ; ylim([-1 6]) ; grid on ; grid minor ; hold on
% 
% uselessplot = plot(0,0); %allows me then to delete c
% pcenter = uselessplot ; c = uselessplot ; v = uselessplot;
% uselesstext = text(0,0,' ');
% t1 = uselesstext ; t2 = uselesstext;
% while 1 %ctrl+c to escape
% center = ginput(1); clc
% delete(pcenter) ; delete(c)
% circonfPos = circle2D(center(1),center(2),r,nbIntervalsForCircle);
% pcenter = plot(center(1),center(2),'r+'); hold on
% c = plot(circonfPos(:,1),circonfPos(:,2),'r-'); hold on
% delete (t2) ; t2 = text(3,5,sentence(flag+1));
% delete(t1)
% t1 = text(3,5.5,['center = ',num2str(center)]); %debug

%% function
function yesorno = do1convex2dPolygon1circle2dIntersect(VerticesNx2,center,r)
%in 2D
%Vertices must be ordoned here in clockwise sense!
%Polygons must be convex!

m = size(VerticesNx2,1);

%Does the circle overlap any of the edges?
%To comment for mode function instead of mode script:
% flag = 0; %First let's assume that no
yesorno = 0;
for i=1:m %vertex 1
%         disp(['i = ',num2str(i)])
    j=i+1; %vertex 2
        if j>m
%             disp('j>m')
            j=1;
        end
%         disp(['j = ',num2str(j)])
    extr = [VerticesNx2(i,:),VerticesNx2(j,:)];
%         disp(['extr = ',num2str(extr)])
%     disp(['extr1 = ',num2str(VerticesNx2(i,:))])
%     disp(['extr2 = ',num2str(VerticesNx2(j,:))])
%     disp(['center = ',num2str(center)])
%     disp(['r = ',num2str(r)])
    result = do1Segment1CircleIntersect(extr,center,r);
%         disp(['result = ',num2str(result)])
    if result
        yesorno = 1; %disp('1Segment1CircleIntersect') %debug
        return
%                 flag = 1; %can be passed as a commentary, just used when commentating
                %the head title function and ending lines, to plot and verify
%                 disp(['circle intersect some edge']) %debug
%                 break
    end ; end

% disp('is called though') %debug

% if flag == 0
if yesorno == 0
%Does the circle strictly lies within the polygon?
%(simpler that nonstrictly lies, thanks to the safety circle area)
%totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html
    outsideFlag = 0; %Let's first assume that yes
    %Construct the orthogonal-to-edge vectors: %TODO : check if they are
    %already stored, so that the construction could occurs only once at the
    %beginning of the simulation
    %+Construct of the starting edge point to the point of interest
    for i=1:m %vertex 1
%         delete(v) ; v = [];
%             disp(['i = ',num2str(i)])
        j=i+1; %vertex 2
        if j>m ; j=1 ; end ; %disp(['j = ',num2str(j)])
        leftOrthoToEdgeVectors(i,:) = [-(VerticesNx2(j,2)-VerticesNx2(i,2)) , VerticesNx2(j,1)-VerticesNx2(i,1)];
        startingEdgeToInterestVectors(i,:) = center-VerticesNx2(i,:);
%         middleEdgePoint = (VerticesNx2(i,:)+VerticesNx2(j,:))/2; %tmp
%         v(end+1) = quiver(middleEdgePoint(1),middleEdgePoint(2),leftOrthoToEdgeVectors(i,1),leftOrthoToEdgeVectors(i,2),'AutoScale','off'); %left vector
%         v(end+1) = quiver(VerticesNx2(i,1),VerticesNx2(i,2),startingEdgeToInterestVectors(i,1),startingEdgeToInterestVectors(i,2),'AutoScale','off'); %left vector
        hold on
        test = dot(leftOrthoToEdgeVectors(i,:),startingEdgeToInterestVectors(i,:));
%         disp(['dot = ',num2str(test)])
%         pause()
        if test > 0 %center lies outside of the ASSUMED CONVEX polygon
            outsideFlag = 1;
            break
        end
    end
    if outsideFlag == 0;
%         flag = 1;
        yesorno = 1; %disp('inside') %debug
        return
%         disp(['inside']) %debug
    %Not intuitive concept -> TODO : find the demo + design an animation

    else
    %     return 0;
%         flag = 0;
        yesorno = 0;
%         disp(['outside']) %debug
    end
end

% end
end

%one iteration of delay might remain in script mode, but normally no delay
%when the function is used as a function and not a script..
