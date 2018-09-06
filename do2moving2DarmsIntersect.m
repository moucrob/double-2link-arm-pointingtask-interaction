% orig1x = 0.57/2 + 0.35/2; %later dont forget to add radius of the arm cylinder
% orig1y = -((0.7-0.345)/2) -0.3;
% orig1 = [orig1x orig1y];
% orig2x = 0.57/2 - 0.35/2; %later dont forget to add radius of the arm cylinder
% orig2y = 0.345 + (0.7-0.345)/2 +0.3;
% orig2 = [orig2x orig2y];
% origx4 = [orig1 orig2]; %x1 y1 x2 y2
% argx8 = [0 0 0 0.345 0.57 0.345 0.57 0]; % m
% rx2 = [1 1]*0.005 ; durx2  = [1 1] ; display = 1;
% dimsArms = [ones(2,1)*0.5 ones(2,1)*0.5];

function yesorno = do2moving2DarmsIntersect(argx8,rx2,dimsArms,orig1,orig2,durx2,display)
%To see whether 2 2D-arms, each coming and backing along a straight line
%between two end points, will intersect or not. Their velocity is such that
%I use a polynomial interpolation (order 3 bc 4 cond, and not higher bc
%velocities through time are low enough) with initial and final velocities = 0.

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

r1 = rx2(1); %radius of the fingerprint (m)
r2 = rx2(2); %m

r = 0.05; %%% TUNABLE %%% radius of the oriented angles (m)
scal = 0.03; %%% TUNABLE %%% length of the two side parts of the arrow heads

dur1 = durx2(1); % total duration including coming and backing moves.
dur2 = durx2(2);

discr = 100; %number of points computed by trajectory
t1 = linspace(0,dur1,discr);
dt1 = t1(2);
t2 = linspace(0,dur2,discr);

a1_1 = dimsArms(1,1); %1st arm
a2_1 = dimsArms(1,2); %1st arm
a1_2 = dimsArms(2,1);
a2_2 = dimsArms(2,2);

posStart1 = [x1start y1start];
posStart2 = [x2start y2start];
posEnd1 = [x1end y1end];
posEnd2 = [x2end y2end];
orig1x = orig1(1) ; orig1y = orig1(2);
orig2x = orig2(1) ; orig2y = orig2(2);

EEs = 0 ; foreArms = 0 ; arms = 0; %flags for bumping
fore1arm2 = 0 ; fore2arm1 = 0;
EE1fore2 = 0 ; EE1arm2 = 0;
EE2fore1 = 0 ; EE2arm1 = 0;

%Trajectories:
% global Q1 Q2 %just for testing something, to remove!
P1 = [] ; P2 = [] ; Q1 = [] ; Q2 = [] ; Elb1 = [] ; Elb2 = [];

%The indented blocks don't work, dunno why...
%     Tbottom = [0 0 1];
%     R10 = eye(2) ; disp(['R10 = ',num2str(R10)])
%     t10 = orig1'; disp(['t10 = ',num2str(t10)])
%     T10 = [[R10,t10];Tbottom]; disp(['T10 = ',num2str(T10)])
%     R20 = -eye(2) ; t20 = orig2'; T20 = [[R20,t20];Tbottom]; disp(['R20 = ',num2str(R20)])
    
if dur1 == dur2 %we can only consider one step of the coming and backing move
    yesorno = 0; %by default
    for i=1:discr
%         disp(['i = ',num2str(i)])
        P1(i,:) = lagrInterp2posCond2nullVelCond([posStart1 posEnd1],dur1,t1(i));
        P2(i,:) = lagrInterp2posCond2nullVelCond([posStart2 posEnd2],dur2,t2(i));
        %Inverse kinematics :
        %but before, frame changement:
        
%             P1 = P1(i,:)' ; P1(end+1) = 1;
%             P2 = P2(i,:)' ; P2(end+1) = 1;
        
        P1fromOrigin1 = P1(i,:) - orig1; %I keep the orientation of the global tablet
        P2fromOrigin2 = P2(i,:) - orig2; %I keep the orientation of the global tablet
        
%             P1fromOrigin1 = T10*P1 ; P1fromOrigin1 = P1fromOrigin1(1:end-1)';
%             P2fromOrigin2 = T20*P2 ; P2fromOrigin2 = P2fromOrigin2(1:end-1)';
        
        Q1(i,:) = IKright2RarmWithoutAlKashi(P1fromOrigin1,a1_1,a2_1); %upper Arm 1, foreArm1
%         disp(['Q1(',num2str(i),',:) = ',num2str(Q1(i,:)*180/pi)])
        Q2(i,:) = IKright2RarmWithoutAlKashi(P2fromOrigin2,a1_2,a2_2);
%         disp(['Q2(',num2str(i),',:) = ',num2str(Q2(i,:)*180/pi)])
        Elb1(i,:) = orig1 + [a1_1*cos(Q1(i,1)) , a1_1*sin(Q1(i,1))]; %Q1(i,1) = q1(t_i)
%         disp(['Elb1(',num2str(i),',:) = ',num2str(Elb1(i,:))])
        Elb2(i,:) = orig2 + [a1_2*cos(Q2(i,1)) , a1_2*sin(Q2(i,1))];
%         disp(['Elb2(',num2str(i),',:) = ',num2str(Elb2(i,:))])
%         posJoints1(i,:) = [posStart1 something posEnd1];
        %A bunch of if bumping...
        if do1Segment1CircleIntersect([Elb2(i,:) P2(i,:)], P1(i,:), r1)
            EE1fore2 = 1; yesorno = 1; break; end
        if do1Segment1CircleIntersect([Elb1(i,:) P1(i,:)], P2(i,:), r2)
            EE2fore1 = 1; yesorno = 1; break; end
        if do1Segment1CircleIntersect([orig2 Elb2(i,:)], P1(i,:), r1)
            EE1arm2 = 1; yesorno = 1; break; end
        if do1Segment1CircleIntersect([orig1 Elb1(i,:)], P2(i,:), r2)
            EE2arm1 = 1; yesorno = 1; break; end
        if do2segmentsIntersect([Elb1(i,:) Elb2(i,:) P1(i,:) P2(i,:)])
            foreArms = 1; yesorno = 1; break; end
        if norm(P2(i,:)-P1(i,:)) < r1+r2
            EEs = 1; yesorno = 1; break; end
        if do2segmentsIntersect([Elb1(i,:) P2(i,:) P1(i,:) Elb2(i,:)])
            fore1arm2 = 1; yesorno = 1; break; end
        if do2segmentsIntersect([P1(i,:) Elb2(i,:) Elb1(i,:) P2(i,:)])
            fore2arm1 = 1; yesorno = 1; break; end %forearm2 upperArm1
        if do2segmentsIntersect([P1(i,:) P2(i,:) Elb1(i,:) Elb2(i,:)])
            arms = 1; yesorno = 1; break; end %upperArms
    end
else %we have to simulate the move from the coming step UNTIL the backing step
    disp('Duration1 different from duration2, case still to be implemented.') %too lazy 
end

%Velocities:
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
    %anchor the plot w.r.t the elbow the most out-of-the-tablet config:
%     Goal1ForMostOut = [0.57 0]; Goal1frame1 = Goal1ForMostOut - orig1;
%     disp(['Goal1frame1 = ',num2str(Goal1frame1)])
%     THET1 = IKright2RarmWithoutAlKashi(Goal1frame1,a1_1,a2_1);
%     disp(['THET1 = ',num2str(THET1*(180/pi))])
%     Goal2ForMostOut = [0 0.345]; Goal2frame2 = Goal2ForMostOut - orig2;
%     disp(['Goal2frame2 = ',num2str(Goal2frame2)])
%     THET2 = IKright2RarmWithoutAlKashi(Goal2frame2,a1_2,a2_2);
%     disp(['THET2 = ',num2str(THET2*(180/pi))])
%     xlim([orig2(1)+a1_2*cos(THET2(1)) orig1(1)+a1_1*cos(THET1(1))])
    xlim([orig2(1)-a1_2 orig1(1)+a1_1])
    zlim([0 inf]) ; ylim([orig1(2) orig2(2)])
%     %debug:
%     X1debug = [0 ; a1_1*cos(THET1(1)) ; a1_1*cos(THET1(1)) + a2_1*cos(sum(THET1))] + orig1(1);
%     Y1debug = [0 ; a1_1*sin(THET1(1)) ; a1_1*sin(THET1(1)) + a2_1*sin(sum(THET1))] + orig1(2);
%     debug1 = plot3(X1debug,Y1debug,zeros(3,1),'k-') ; hold on
%     X2debug = [0 ; a1_2*cos(THET2(1)) ; a1_2*cos(THET2(1)) + a2_2*cos(sum(THET2))] + orig2(1);
%     Y2debug = [0 ; a1_2*sin(THET2(1)) ; a1_2*sin(THET2(1)) + a2_2*sin(sum(THET2))] + orig2(2);
%     debug2 = plot3(X2debug,Y2debug,zeros(3,1),'k-') ; hold on
    
    %paths:
    quiver3(x1start,y1start,0,(x1end-x1start),(y1end-y1start),0,'--g','AutoScale','off')
    grid on ; grid minor ; xlabel('x (cm)') ; ylabel('y (cm)'), zlabel ('vel (m/s)') ; hold on
    quiver3(x2start,y2start,0,(x2end-x2start),(y2end-y2start),0,'--b','AutoScale','off') ; hold on
    
    %interactive plot:
    nbIntervalsForCircle = 50 ; nbPointsForCircle = nbIntervalsForCircle + 1;
    uselessplot = plot3(0,0,0) ; circle1 = uselessplot; circle2 = uselessplot;
    segment1 = uselessplot ; segment2 = uselessplot;
    Zc = zeros(nbPointsForCircle,1); obs = 1; %time of pause
    Z2 = zeros(2,1);
    arm1 = uselessplot ; fore1 = uselessplot;
    arm2 = uselessplot ; fore2 = uselessplot;
    angle1_p1 = uselessplot ; angle2_p1 = uselessplot;
    angle1_p2 = uselessplot ; angle2_p2 = uselessplot;
    t = text(0,0,0,'useless');
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
        delete(circle1) ; circle1 = plot3(tmp1(:,1),tmp1(:,2),Zc,'g'); hold on
        delete(circle2) ; circle2 = plot3(tmp2(:,1),tmp2(:,2),Zc,'b'); hold on
        delete(arm1) ; arm1 = plot3([orig1x;Elb1(j,1)],[orig1y;Elb1(j,2)],Z2,'g'); hold on
            %debug:
%             disp(['length arm1 = ',num2str(norm(Elb1(j,:)-orig1))])
        delete(arm2) ; arm2 = plot3([orig2x;Elb2(j,1)],[orig2y;Elb2(j,2)],Z2,'b'); hold on
            %debug:
%             disp(['length arm2 = ',num2str(norm(Elb2(j,:)-orig2))])
        delete(fore1) ; fore1 = plot3([Elb1(j,1);P1(j,1)],[Elb1(j,2);P1(j,2)],Z2,'g'); hold on
%         fore1 = plot3([Elb1(j,1);Elb1(j,1)+cos(sum(Q1(j,:)))],[Elb1(j,2);Elb1(j,2)+sin(sum(Q1(j,:)))],Z2,'g'); hold on    
            %debug: this shows that length of the forearms vary according
            %to the time !!!!!! :(
%             disp(['length fore1 = ',num2str(norm(Elb1(j,:)-P1(j,:)))])
        delete(fore2) ; fore2 = plot3([Elb2(j,1);P2(j,1)],[Elb2(j,2);P2(j,2)],Z2,'b'); hold on
%         fore2 = plot3([Elb2(j,1);Elb2(j,1)+cos(sum(Q2(j,:)))],[Elb2(j,2);Elb2(j,2)+sin(sum(Q2(j,:)))],Z2,'b'); hold on      
            %debug:
%             disp(['length fore2 = ',num2str(norm(Elb2(j,:)-P2(j,:)))])
            
        delete(angle1_p1);
        angle1_p1 = plot_arc(0,Q1(j,1),orig1,r,0,'g',scal);
        delete(angle2_p1);
        angle2_p1 = plot_arc(Q1(j,1),sum(Q1(j,:)),Elb1(j,:),r,0,'g',scal);
        delete(angle1_p2);
        angle1_p2 = plot_arc(0,Q2(j,1),orig2,r,0,'b',scal);
        delete(angle2_p2);
        angle2_p2 = plot_arc(Q2(j,1),sum(Q2(j,:)),Elb2(j,:),r,0,'b',scal);
        
        delete(t)
        t = text(ones(4,1)*0.6,[0.65 0.55 0.45 0.35],zeros(4,1), ...
                 {['q_{1_{player_{1}}}^{R_{0}} = ',num2str(Q1(j,1)*180/pi),'°'], ...
                  ['q_{2_{player_{1}}}^{R_{0}} = ',num2str(Q1(j,2)*180/pi),'°'], ...
                  ['q_{1_{player_{2}}}^{R_{0}} = ',num2str(Q2(j,1)*180/pi),'°'], ...
                  ['q_{2_{player_{2}}}^{R_{0}} = ',num2str(Q2(j,2)*180/pi),'°'], ...
                 }, ...
                'FontSize',8);
        
        drawnow
%         pause(1)
%         if j==1 ; pause(20) ; end %debug
            if (j==i) ... %we arrive at the bumping moment
            && (yesorno == 1)
                if EE1fore2;set(circle1,'Color','r'); set(fore2,'Color','r'); end
                if EE2fore1;set(circle2,'Color','r'); set(fore1,'Color','r'); end
                if EE1arm2;set(circle1,'Color','r'); set(arm2,'Color','r'); end
                if EE2arm1;set(circle2,'Color','r'); set(arm1,'Color','r'); end
                if foreArms;set(fore1,'Color','r'); set(fore2,'Color','r'); end
                if EEs; set(circle1,'Color','r'); set(circle2,'Color','r'); end
                if fore1arm2;set(fore1,'Color','r'); set(arm2,'Color','r'); end
                if fore2arm1;set(fore2,'Color','r'); set(arm1,'Color','r'); end
                if arms;set(arm1,'Color','r'); set(arm2,'Color','r'); end
                pause(5)
            end
    end

    plotTime = toc;
    disp(['The plot took ',num2str(plotTime-obs),' s.'])
end

end