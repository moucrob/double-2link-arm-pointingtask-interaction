%To use this script as a script instead of a func to see the workspace:
r = 1;
nbIntervalsForCircle = 50;
set(0,'DefaultFigureWindowStyle','docked')
figure;
set(0,'DefaultFigureWindowStyle','normal')

uselessplot = plot(0,0); %allows me then to delete c
point1 = uselessplot ; point2 = uselessplot ; seg = uselessplot;
pcenter = uselessplot ; c = uselessplot ; points = uselessplot;
uselesstext = text(0,0,' '); t = uselesstext ;
axis square ; xlim([0 6]) ; ylim([0 6]) ; grid on ; grid minor ; hold on

while 1
argx4 = [];
delete(t) ; t = text(2,5, 'Click on the starting point of the segment');
argx4(1,1:2) = ginput(1);
delete(point1) ; point1 = plot(argx4(1),argx4(2),'s','MarkerSize',10);
delete(point2) ; delete(seg);
delete(pcenter) ; delete(c) ; delete(points);

delete(t) ; t = text(2,5, 'Click on the ending point of the segment');
argx4(1,3:4) = ginput(1);
point2 = plot(argx4(3),argx4(4),'s','MarkerSize',10);
seg = plot([argx4(1);argx4(3)],[argx4(2);argx4(4)]) ; hold on

delete(t) ; t = text(2,5, 'Click on the center of the circle');
Center = ginput(1);
pcenter = plot(Center(1),Center(2),'r+') ; hold on 
circonfPos = circle2D(Center(1),Center(2),r,nbIntervalsForCircle);
c = plot(circonfPos(:,1),circonfPos(:,2),'r-') ; hold on

%https://math.stackexchange.com/questions/2837/how-to-tell-if-a-line-segment-intersects-with-a-circle
% function yesorno = do1Segment1CircleIntersect(argx4,Center,r)
%input
xstart = argx4(1);
ystart = argx4(2);
xend = argx4(3);
yend = argx4(4);
xC = Center(1) ; yC = Center(2);

%formatting:
P = [xstart ; ystart];
Q = [xend ; yend];
C = [xC ; yC];

%I assume the line is much longer than the diameter of the circle,
%otherwise the problem is not completely explored !!
sentence = strcat("This function is designed such that \n the circle's ", ...
                  "diameter is lower than the segment length. \n Here", ...
                  "this is not the case, \n so one should have a look in ", ...
                  "the snippet and resolve the problem in a complete way.");
length = norm(P-Q);
if norm(P-Q) < 2*r
    clc
    fprintf(sentence)
    disp(' ')
    disp(['(length = ',num2str(length),')'])
end

%Mathematica-given solution:
%(each point of a line passing by P and Q is described by one unique
%parameter t s.t tP+(1-t)Q belongs to the line, while if t is borned
%between 0 and 1, then tP+(1-t)Q belongs to the segment between P and Q.)
underSqrt = P(1)^2*(-Q(2)^2) - 2*P(1)*Q(1) + 2*P(2)*P(1)*Q(1)*Q(2) ...
          - P(2)^2*Q(1)^2 - 2*P(2)*Q(2) + P(1)^2 + P(2)^2 + Q(1)^2 + Q(2)^2;
if underSqrt < 0
    yesorno = 0;
else
    num1 = -P'*Q;
    num2 = Q'*Q;
    denom = 2*num1 + P'*P + num2;
    t1 = (num1 - sqrt(underSqrt) + num2)/denom; %first segment-parameter root
    t2 = (num1 + sqrt(underSqrt) + num2)/denom;
    if ((0 <= t1) && (t1 <= 1)) || ((0 <= t2) && (t2 <= 1)); yesorno = 1;
    else; yesorno = 0; end
end

%To verify this function isolated of any bigger main script, uncomment:
R1 = t1*P + (1-t1)*Q ; R2 = t2*P + (1-t2)*Q;
points = plot([R1(1) R2(1)],[R1(2) R2(2)],'.','MarkerSize',20);
title(['yesorno = ', num2str(yesorno)])

% end

end