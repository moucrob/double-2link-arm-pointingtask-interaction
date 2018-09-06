%% Workspace for verifying the function (comment the function title and ending lines):
% r = 1;
% nbIntervalsForCircle = 50;
% set(0,'DefaultFigureWindowStyle','docked')
% figure;
% set(0,'DefaultFigureWindowStyle','normal')
% 
% uselessplot = plot(0,0); %allows me then to delete c
% point1 = uselessplot ; point2 = uselessplot ; seg = uselessplot;
% pcenter = uselessplot ; c = uselessplot ; points = uselessplot;
% uselesstext = text(0,0,' '); t = uselesstext ;
% axis square ; xlim([0 6]) ; ylim([0 6]) ; grid on ; grid minor ; hold on
% 
% while 1
% argx4 = [];
% delete(t) ; t = text(2,5, 'Click on the starting point of the segment');
% argx4(1,1:2) = ginput(1);
% delete(point1) ; point1 = plot(argx4(1),argx4(2),'s','MarkerSize',10);
% delete(point2) ; delete(seg);
% delete(pcenter) ; delete(c) ; delete(points);
% 
% delete(t) ; t = text(2,5, 'Click on the ending point of the segment');
% argx4(1,3:4) = ginput(1);
% point2 = plot(argx4(3),argx4(4),'s','MarkerSize',10);
% seg = plot([argx4(1);argx4(3)],[argx4(2);argx4(4)]) ; hold on
% 
% delete(t) ; t = text(2,5, 'Click on the center of the circle');
% Center = ginput(1);
% pcenter = plot(Center(1),Center(2),'r+') ; hold on 
% circonfPos = circle2D(Center(1),Center(2),r,nbIntervalsForCircle);
% c = plot(circonfPos(:,1),circonfPos(:,2),'r-') ; hold on

%% Function 
function yesorno = do1Segment1CircleIntersect(argx4,Center,r)
%(include also the case where circle contains the segment)
%input syntax :
%xstart = argx4(1);
%ystart = argx4(2);
%xend = argx4(3);
%yend = argx4(4);
%xC = Center(1) ; yC = Center(2);

%     disp(['distCtoSeg = ',num2str( distPointSegment(argx4,Center) )])
%     disp(['radius',num2str(r)])
if distPointSegment(argx4,Center) <= r ; yesorno = 1 ; else ; yesorno = 0 ; end

%     title(['yesorno = ', num2str(yesorno)])

end

% end