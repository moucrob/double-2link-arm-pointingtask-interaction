% argx8 = [0 0 0 1 1 1 1 0]*30 % cm
% rx2 = [1 1]*0.5 ; durx2  = [1 1] ; display = 1;

function yesorno = do2movingCirclesIntersect(argx8,rx2,durx2,display)
%To see whether 2 circles, each coming and backing along a straight line between two
%end points, will intersect or not. Their velocity is such that I use a
%polynomial interpolation with initial and final velocities = 0.

%We assume that reaction time of the 2 players are same hence they leave
%the standby / piece location positions at exactly the same time.

clf ; clc

x1start = argx8(1);
y1start = argx8(2);
x1end = argx8(5);
y1end = argx8(6);
x2start = argx8(3);
y2start = argx8(4);
x2end = argx8(7);
y2end = argx8(8);

r1 = rx2(1); %cm
r2 = rx2(2); %cm

dur1 = durx2(1); % total duration including coming and backing moves.
dur2 = durx2(2);

discr = 100; %number of points computed by trajectory
t1 = linspace(0,dur1,discr);
dt1 = t1(2);
t2 = linspace(0,dur2,discr);

%Trajectories:
P1 = [] ; P2 = [];
if dur1 == dur2 %we can only consider one step of the coming and backing move
    yesorno = 0; %by default
    for i=1:discr
        P1(i,:) = lagrInterp2posCond2nullVelCond([x1start y1start x1end y1end],dur1,t1(i));
        P2(i,:) = lagrInterp2posCond2nullVelCond([x2start y2start x2end y2end],dur2,t2(i));
        if norm(P2(i,:)-P1(i,:)) < r1+r2; yesorno = 1; break; end %or return
    end
else %we have to simulate the move from the coming step UNTIL the backing step
    disp('Duration1 different from duration2, case still to be implemented.') %too lazy 
end

%Velocities
V1 = [] ; V2 = [];
V1(:,1) = centered2ndOrderDerivativeNullBoundaries(P1(:,1),dt1);
V1(:,2) = centered2ndOrderDerivativeNullBoundaries(P1(:,2),dt1);
V2(:,1) = centered2ndOrderDerivativeNullBoundaries(P2(:,1),dt1);
V2(:,2) = centered2ndOrderDerivativeNullBoundaries(P2(:,2),dt1);
normV1 = sqrt(sum(V1.^2,2)); normV2 = sqrt(sum(V2.^2,2));

%% plots
if display
    tic 

    %tablet boundaries:
    plot3([0;0.57;0.57;0;0],[0;0;0.345;0.345;0],zeros(5,1),'k-','LineWidth',1) ; hold on
    
    % plot3(P1(:,1), P1(:,2),zeros(discr,1)) ; grid on ; grid minor
    %Ok this demonstrates those traj are straight line.
    %We can then allow ourselves to only pot a straigt arrow (it fasteners)
    quiver3(x1start,y1start,0,(x1end-x1start),(y1end-y1start),0,'--g','AutoScale','off')
    grid on ; grid minor ; xlabel('x (cm)') ; ylabel('y (cm)'), zlabel ('vel (m/s)') ; hold on

    % axis([min(argx8(1),argx8(5)) max(argx8(3),argx8(7)) ...
    %       min(argx8(2),argx8(6)) max(argx8(4),argx8(8))])
    ylim([0 0.345])

    % plot3(P2(:,1), P2(:,2),zeros(discr,1))
    quiver3(x2start,y2start,0,(x2end-x2start),(y2end-y2start),0,'--b','AutoScale','off') ; hold on

    %interactive plot:
    nbIntervalsForCircle = 50 ; nbPointsForCircle = nbIntervalsForCircle + 1;
    uselessplot = plot3(0,0,0) ; circle1 = uselessplot; circle2 = uselessplot;
    segment1 = uselessplot ; segment2 = uselessplot;
    Z = zeros(nbPointsForCircle,1); obs = 1; %time of pause
    for j=1:i
        if j>=2
            curvePortion1 = [[P1(j,1) ; P1(j-1,1)] [P1(j,2) ; P1(j-1,2)] [normV1(j) ; normV1(j-1)]];
            curvePortion2 = [[P2(j,1) ; P2(j-1,1)] [P2(j,2) ; P2(j-1,2)] [normV2(j) ; normV2(j-1)]];
            vert1 = [[P1(j,1) ; P1(j,1)] [P1(j,2) ; P1(j,2)] [0 ; normV1(j)]];
            vert2 = [[P2(j,1) ; P2(j,1)] [P2(j,2) ; P2(j,2)] [0 ; normV2(j)]];
            delete(segment1)
            segment1 = plot3(vert1(:,1),vert1(:,2),vert1(:,3),':g'); hold on
            delete(segment2)
            segment2 = plot3(vert2(:,1),vert2(:,2),vert2(:,3),':b'); hold on
            plot3(curvePortion1(:,1),curvePortion1(:,2),curvePortion1(:,3),'-g') ; hold on
            plot3(curvePortion2(:,1),curvePortion2(:,2),curvePortion2(:,3),'-b') ; hold on
        end
        tmp1 = circle2D(P1(j,1),P1(j,2),r1,nbIntervalsForCircle);
        tmp2 = circle2D(P2(j,1),P2(j,2),r2,nbIntervalsForCircle);
        delete(circle1) ; circle1 = plot3(tmp1(:,1),tmp1(:,2),Z,'g'); hold on
        delete(circle2) ; circle2 = plot3(tmp2(:,1),tmp2(:,2),Z,'b'); hold on
        drawnow
            if (j==i) ... %we arrive at the bumping moment
            && (yesorno == 1)
                set(circle1,'Color','r') ; set(circle2,'Color','r') ; pause(obs)
            end
    end

    plotTime = toc;
    disp(['The plot took ',num2str(plotTime-obs),' s.'])
end

end