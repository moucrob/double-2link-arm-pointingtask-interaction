% clear all  %V1{i,1} assignment might require this clear... 
% 
% orig1x = 0.57/2 + 0.35/2; %later dont forget to add radius of the arm cylinder
% orig1y = -((0.7-0.345)/2) -0.3;
% orig1 = [orig1x orig1y];
% orig2x = 0.57/2 - 0.35/2; %later dont forget to add radius of the arm cylinder
% orig2y = 0.345 + (0.7-0.345)/2 +0.3;
% orig2 = [orig2x orig2y];
% origx4 = [orig1 orig2]; %x1 y1 x2 y2
% argx8 = [0.57+0.02 0.345/2 0 0.345/2 -0.02 0.345/2 0.57 0]; %(m) (logic below)
% rx2 = [1 1]*0.005 ; durx2  = [1 1] ; display = 1;
% 
% %VISION ADDITION:
% origVision1x = 0.57/2 ; origVision2x = origVision1x;
% origVision1z = 0.40 ; origVision2z = origVision1z;
% origVision1y = orig1y + 0.094;
% origVision2y = orig2y - 0.094;
% origVision1 = [origVision1x origVision1y origVision1z];
% origVision2 = [origVision2x origVision2y origVision2z];
% 
% dimsArms = [ones(2,1)*0.32 ...
%             ones(2,1)*0.4 ones(2,1)*0.075 ones(2,1)*0.06 ...
%             ones(2,1)*0.1 ones(2,2)*2*rx2(1)*cos(pi/4) ]; %(m) (logic below)
% [whichAngle,lenFore,projectedLenFing] = whichDesignAngleRotY(dimsArms(1,4),dimsArms(1,5),dimsArms(1,7),0.5,0);
% dimsArms(:,2) = ones(2,1)*lenFore;
% [epsilonTot,O2toEEforIK]= perpendicularDecenteringEE(dimsArms(1,3)/2,lenFore+projectedLenFing);
% [epsilonFore,O2toMETforIK]= perpendicularDecenteringEE(dimsArms(1,3)/2,lenFore);

function yesorno = do2_2Dmoving_3DforeArms_decentered3DfingersIntersect(argx8,rx2,dimsArms,orig1,orig2,origVision1,origVision2,durx2,whichAngle,lenFore,projectedLenFing,epsilonTot,epsilonFore,O2toEEforIK,O2toMETforIK,display)
%To see whether 2 2D-arms, each coming and backing along a straight line
%between two end points, will intersect and/or occlude or not. Their velocity is such that
%I use a polynomial interpolation (order 3 bc 4 cond, and not higher bc
%velocities through time are low enough) with initial and final velocities = 0.

%We assume that reaction time of the 2 players are same hence they leave
%the standby / piece location positions at exactly the same time.

clf ; %clc

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

a1_1 = dimsArms(1,1); %length 1st arm
a2_1 = dimsArms(1,2); %length 1st forearm
%In the plane of the table:
b2_1 = dimsArms(1,3); %width of the squared base of the cylinder (1st arm)
c2_1 = dimsArms(1,4); %height of the squared base of the cylinder (1st arm)
a3_1 = dimsArms(1,5); %length 1st finger
b3_1 = dimsArms(1,6); %width 
c3_1 = dimsArms(1,7); %height

a1_2 = dimsArms(2,1);
a2_2 = dimsArms(2,2) ; b2_2 = dimsArms(2,3); c2_2 = dimsArms(2,4);
a3_2 = dimsArms(2,5) ; b3_2 = dimsArms(2,6) ; c3_2 = dimsArms(2,7);

alti = c2_1/2; %altitude of the bones

templateForeX = repmat(a2_1*[0 ; 1 ; 1 ; 0],2,1); %repeat the vector. Putting ,2,1 or ,1,2
%is quite unpredictable, didnt understood the logic, so I prefer not comment this way...
templateForeY = repmat((b2_1/2)*[-1 ; -1 ; 1 ; 1],2,1);
templateForeZ = c2_1*[zeros(4,1);ones(4,1)];
templateForeV = [templateForeX,templateForeY,templateForeZ];
ForeToBeTransformed = [templateForeV';ones(1,size(templateForeV,1))];
templateForeF = [1 2 6 5 ; ... %Face1 connects vertices (rows of V1) number 1,2,5,6
                 2 3 7 6 ; ... 
                 3 4 8 7 ; ...
                 4 1 5 8 ; ...
                 1 2 3 4 ; ...
                 5 6 7 8 ];
Tz = @(theta,t3x1) ...
    [cos(theta) -sin(theta) 0 t3x1(1);
     sin(theta) cos(theta)  0 t3x1(2);
          0          0      1 t3x1(3);
          0          0      0    1   ];
      
templateFingX = repmat(a3_1*[0 ; 1 ; 1 ; 0],2,1); %repeat the vector. Putting ,2,1 or ,1,2
%is quite unpredictable, didnt understood the logic, so I prefer not comment this way...
templateFingY = repmat((b3_1/2)*[-1 ; -1 ; 1 ; 1],2,1);
templateFingZ = c3_1*[zeros(4,1);ones(4,1)];
templateFingV = [templateFingX,templateFingY,templateFingZ];
FingToBeTransformed = [templateFingV';ones(1,size(templateFingV,1))];
templateFingF = templateForeF;
Ry = @(phi) ...
    [cos(phi) 0 sin(phi); ...
        0     1    0    ; ...
    -sin(phi) 0 cos(phi)];
templateFingV = Ry(whichAngle)*templateFingV' ; templateFingV = templateFingV';
FingToBeTransformed = [templateFingV' ; ones(1,size(templateFingV,1))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TO DO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterationsAllowed = 2; %IMPORTANT PARAMETER TO STUDY 

posStart1 = [x1start y1start] ; posStart2 = [x2start y2start];
posEnd1 = [x1end y1end] ; posEnd2 = [x2end y2end];
orig1x = orig1(1) ; orig1y = orig1(2);
orig2x = orig2(1) ; orig2y = orig2(2);

EEs = 0 ; foreArms = 0 ; arms = 0; %flags for bumping
fore1arm2 = 0 ; fore2arm1 = 0;
EE1fore2 = 0 ; EE1arm2 = 0;
EE2fore1 = 0 ; EE2arm1 = 0;
%additions:
fingers = 0;
fing1arm2 = 0 ; fing2arm1 = 0;
fing1fore2 = 0 ; fing2fore1 = 0;
% EE1fing2 = 0 ; EE2fing1 = 0;
%again additions:
vision1fing2 = 0 ; vision2fing1 = 0;
vision1fore2 = 0 ; vision2fore1 = 0;

%Trajectories:
% global Q1 Q2 %just for testing something, to remove!
P1 = [] ; P2 = [];
Q1 = [] ; Q2 = [];
Elb1 = [] ; Elb2 = [];
V1 = [] ; V2 = [];
Met1 = [] ; Met2 = []; %metacarp
V3 = [] ; V4 = [];

%The indented blocks don't work, dunno why...
%     Tbottom = [0 0 1];
%     R10 = eye(2) ; disp(['R10 = ',num2str(R10)])
%     t10 = orig1'; disp(['t10 = ',num2str(t10)])
%     T10 = [[R10,t10];Tbottom]; disp(['T10 = ',num2str(T10)])
%     R20 = -eye(2) ; t20 = orig2'; T20 = [[R20,t20];Tbottom]; disp(['R20 = ',num2str(R20)])
    
if dur1 == dur2 %we can only consider one step of the coming and backing move
    yesorno = 0 ; collisionOrNo = 0 ; occlusionOrNo = 0; %by default
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
        
%         Q1(i,:) = IKright2RarmWithoutAlKashi(P1fromOrigin1,a1_1,a2_1+projectedLenFing); %upper Arm 1, foreArm1
        Q1(i,:) = IKright2RarmWithoutAlKashi(P1fromOrigin1,a1_1,O2toEEforIK);
        Q1(i,2) = Q1(i,2)-epsilonTot;
%         disp(['Q1(',num2str(i),',:) = ',num2str(Q1(i,:)*180/pi)])
%         Q2(i,:) = IKright2RarmWithoutAlKashi(P2fromOrigin2,a1_2,a2_2+projectedLenFing);
        Q2(i,:) = IKright2RarmWithoutAlKashi(P2fromOrigin2,a1_2,O2toEEforIK);
        Q2(i,2) = Q2(i,2)-epsilonTot;
%         disp(['Q2(',num2str(i),',:) = ',num2str(Q2(i,:)*180/pi)])
        Elb1(i,:) = orig1 + [a1_1*cos(Q1(i,1)) , a1_1*sin(Q1(i,1))]; %Q1(i,1) = q1(t_i)
%         disp(['Elb1(',num2str(i),',:) = ',num2str(Elb1(i,:))])
        Elb2(i,:) = orig2 + [a1_2*cos(Q2(i,1)) , a1_2*sin(Q2(i,1))];
%         disp(['Elb2(',num2str(i),',:) = ',num2str(Elb2(i,:))])
%         posJoints1(i,:) = [posStart1 something posEnd1];
        Met1(i,:) = Elb1(i,:) + [O2toMETforIK*cos(sum(Q1(i,:))+epsilonFore) , O2toMETforIK*sin(sum(Q1(i,:))+epsilonFore)];
        Met2(i,:) = Elb2(i,:) + [O2toMETforIK*cos(sum(Q2(i,:))+epsilonFore) , O2toMETforIK*sin(sum(Q2(i,:))+epsilonFore)];
        
        %Now I have to "represent" my foreArm cubes by using t3x1 as Elb_1,2 (2x1 only):
        toTrim1 = Tz(sum(Q1(i,:)),[Elb1(i,:),0])*ForeToBeTransformed;
        V1{i,1} = toTrim1(1:3,:)';
        toTrim2 = Tz(sum(Q2(i,:)),[Elb2(i,:),0])*ForeToBeTransformed;
        V2{i,1} = toTrim2(1:3,:)';
        Fore1obj = patch('Vertices',V1{i,1},'Faces',templateForeF,'visible','off');
        Fore2obj = patch('Vertices',V2{i,1},'Faces',templateForeF,'visible','off');
        %additions:
        toTrim3 = Tz(sum(Q1(i,:)),[Met1(i,:),c2_1-c3_1/2])*FingToBeTransformed;
        V3{i,1} = toTrim3(1:3,:)';
        toTrim4 = Tz(sum(Q2(i,:)),[Met2(i,:),c2_2-c3_2/2])*FingToBeTransformed;
        V4{i,1} = toTrim4(1:3,:)';
        Fing1obj = patch('Vertices',V3{i,1},'Faces',templateForeF,'visible','off');
        Fing2obj = patch('Vertices',V4{i,1},'Faces',templateForeF,'visible','off');

        %A bunch of if bumping... try to hierarchize by most probable to
        %minimize loss of time
        while 1 %try all the possible collisions
            if do1convex2dPolygon1circle2dIntersect(V2{i,1}(1:4,1:2), P1(i,:), r1)
                EE1fore2 = 1; collisionOrNo = 10; break; end
            if do1convex2dPolygon1circle2dIntersect(V1{i,1}(1:4,1:2), P2(i,:), r2)
                EE2fore1 = 1; collisionOrNo = 10; break; end
            if GJK(Fing1obj,Fing2obj,iterationsAllowed)
                fingers = 1; collisionOrNo = 10; break; end
            if GJK(Fing1obj,Fore2obj,iterationsAllowed)
                fing1fore2 = 1; collisionOrNo = 10; break; end
            if GJK(Fing2obj,Fore1obj,iterationsAllowed)
                fing2fore1 = 1; collisionOrNo = 10; break; end
            if GJK(Fore1obj,Fore2obj,iterationsAllowed)
                foreArms = 1; collisionOrNo = 10; break; end
    %         if do1Segment1CircleIntersect([orig2 Elb2(i,:)], P1(i,:), r1)
    %             EE1arm2 = 1; collisionOrNo = 1; break; end
    %         if do1Segment1CircleIntersect([orig1 Elb1(i,:)], P2(i,:), r2)
    %             EE2arm1 = 1; collisionOrNo = 1; break; end
            if norm(P2(i,:)-P1(i,:)) < r1+r2
                EEs = 1; collisionOrNo = 10; break; end
            if do1segment3d1ConvexPolyhedronIntersect([orig2 alti ; Elb2(i,:) alti] , V3{i,1},templateForeF,0)
                fing1arm2 = 1; collisionOrNo = 10; break; end
            if do1segment3d1ConvexPolyhedronIntersect([orig1 alti ; Elb1(i,:) alti] , V4{i,1},templateForeF,0)
                fing2arm1 = 1; collisionOrNo = 10; break; end
            if do1segment3d1ConvexPolyhedronIntersect([orig2 alti ; Elb2(i,:) alti] , V1{i,1},templateForeF,0)
                fore1arm2 = 1; collisionOrNo = 10; break; end
            if do1segment3d1ConvexPolyhedronIntersect([orig1 alti ; Elb1(i,:) alti] , V2{i,1},templateForeF,0)
                fore2arm1 = 1; collisionOrNo = 10; break; end %forearm2 upperArm1
            if do2segmentsIntersect([orig1 orig2 Elb1(i,:) Elb2(i,:)])
                arms = 1; collisionOrNo = 10; break; end %upperArms
            break %if no collision happened
        end
        while 1 %try all the possible occlusions
            if do1segment3d1ConvexPolyhedronIntersect([origVision1 ; P1(i,:) 0] , V4{i,1},templateForeF,0) %I assume the players do not occlude themselves, hence choosing appropriated movement, which is completetely wrong since I choose for them the most obvious movement!...
                vision1fing2 = 1; occlusionOrNo = 21; break; end % ...1 means this is player1 who is occluded, 2... means some player is occluded
            if do1segment3d1ConvexPolyhedronIntersect([origVision2 ; P2(i,:) 0] , V3{i,1},templateForeF,0)
                vision2fing1 = 1; occlusionOrNo = 22; break; end
            if do1segment3d1ConvexPolyhedronIntersect([origVision1 ; P1(i,:) 0] , V2{i,1},templateForeF,0)
                vision1fore2 = 1; occlusionOrNo = 21; break; end
            if do1segment3d1ConvexPolyhedronIntersect([origVision2 ; P2(i,:) 0] , V1{i,1},templateForeF,0)
                vision2fore1 = 1; occlusionOrNo = 22; break; end
            break %if no occlusion happened
        end
    end
else %we have to simulate the move from the coming step UNTIL the backing step
    disp('Duration1 different from duration2, case still to be implemented.') %too lazy 
end

%Velocities:
v1 = [] ; v2 = [];
v1(:,1) = centered2ndOrderDerivativeNullBoundaries(P1(:,1),dt1);
v1(:,2) = centered2ndOrderDerivativeNullBoundaries(P1(:,2),dt1);
v2(:,1) = centered2ndOrderDerivativeNullBoundaries(P2(:,1),dt1);
v2(:,2) = centered2ndOrderDerivativeNullBoundaries(P2(:,2),dt1);
normv1 = sqrt(sum(v1.^2,2)); normv2 = sqrt(sum(v2.^2,2));

%% plots
if display
    tic

    mult = 8000; %otherwise it's too much complicated to plot in 3D on two
    %differently scaled ZAxis...
    
    txtftsize = 10; %text fontsize
    
    %transparency:
    alpha = 0.5;
    
    %tablet boundaries:
    plot3([0;0.57;0.57;0;0],[0;0;0.345;0.345;0],zeros(5,1),'k-','LineWidth',1) ; hold on
%     set(0,'DefaultFigureWindowStyle','docked')
%     set(0,'DefaultFigureWindowStyle','normal')
    set(gca,'DataAspectRatio',[1 1 1])
    view(3) ; camorbit(20,0,'data',[0 0 1])
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

    grid on ; grid minor ; xlabel('x (m)') ; ylabel('y (m)'), zlabel ([num2str(mult),' \times vel (m/s)']) ; hold on
    
%         ax1 = gca;
%         ax2 = ax1 ; set(ax2,'Visible','off')

    %paths:
    quiver3(x1start,y1start,0,(x1end-x1start),(y1end-y1start),0,'--g','AutoScale','off')
    quiver3(x2start,y2start,0,(x2end-x2start),(y2end-y2start),0,'--b','AutoScale','off') ; hold on
    
    %interactive plot:
    nbIntervalsForCircle = 50 ; nbPointsForCircle = nbIntervalsForCircle + 1;
    uselessplot = plot3(0,0,0) ; circle1 = uselessplot; circle2 = uselessplot;
    segment1 = uselessplot ; segment2 = uselessplot;
    Zc = zeros(nbPointsForCircle,1); obs = 1; %time of pause
    Z2 = ones(2,1)*(c2_1/2);
    arm1 = uselessplot ; fore1 = uselessplot;
    arm2 = uselessplot ; fore2 = uselessplot;
    angle1_p1 = uselessplot ; angle2_p1 = uselessplot;
    angle1_p2 = uselessplot ; angle2_p2 = uselessplot;
    finger1 = uselessplot ; finger2 = uselessplot;
    vision1 = uselessplot ; vision2 = uselessplot;
    t1 = text(0,0,0,' ') ; t2 = text(0,0,0,' ');
    for j=1:i
        if j>=2
            curvePortion1 = [[P1(j,1) ; P1(j-1,1)] [P1(j,2) ; P1(j-1,2)] mult*[normv1(j) ; normv1(j-1)]];
            curvePortion2 = [[P2(j,1) ; P2(j-1,1)] [P2(j,2) ; P2(j-1,2)] mult*[normv2(j) ; normv2(j-1)]];
            vert1 = [[P1(j,1) ; P1(j,1)] [P1(j,2) ; P1(j,2)] mult*[0 ; normv1(j)]];
            vert2 = [[P2(j,1) ; P2(j,1)] [P2(j,2) ; P2(j,2)] mult*[0 ; normv2(j)]];
            delete(segment1)
            segment1 = plot3(vert1(:,1),vert1(:,2),vert1(:,3),':g'); hold on %do not misunderstand v = velocity, V = vertices, vert = vertical segment
            delete(segment2)
            segment2 = plot3(vert2(:,1),vert2(:,2),vert2(:,3),':b'); hold on
            plot3(curvePortion1(:,1),curvePortion1(:,2),curvePortion1(:,3),'-g') ; hold on
            plot3(curvePortion2(:,1),curvePortion2(:,2),curvePortion2(:,3),'-b') ; hold on
        end
        tmp1 = circle2D(P1(j,1),P1(j,2),r1,nbIntervalsForCircle);
        tmp2 = circle2D(P2(j,1),P2(j,2),r2,nbIntervalsForCircle);
        delete(circle1) ; circle1 = plot3(tmp1(:,1),tmp1(:,2),Zc,'g'); hold on
        delete(circle2) ; circle2 = plot3(tmp2(:,1),tmp2(:,2),Zc,'b'); hold on
        delete(vision1) ; vision1 = plot3([origVision1(1) P1(j,1)],[origVision1(2) P1(j,2)],[origVision1(3) 0],'g-.'); hold on
        delete(vision2) ; vision2 = plot3([origVision2(1) P2(j,1)],[origVision2(2) P2(j,2)],[origVision2(3) 0],'b-.'); hold on
        delete(arm1) ; arm1 = plot3([orig1x;Elb1(j,1)],[orig1y;Elb1(j,2)],Z2,'g'); hold on
            %debug:
%             disp(['length arm1 = ',num2str(norm(Elb1(j,:)-orig1))])
        delete(arm2) ; arm2 = plot3([orig2x;Elb2(j,1)],[orig2y;Elb2(j,2)],Z2,'b'); hold on
            %debug:
%             disp(['length arm2 = ',num2str(norm(Elb2(j,:)-orig2))])
        delete(fore1) ; fore1 = patch('Vertices',V1{j,1},'Faces',templateForeF,'FaceColor','g','FaceAlpha',alpha); hold on
%         fore1 = plot3([Elb1(j,1);Elb1(j,1)+cos(sum(Q1(j,:)))],[Elb1(j,2);Elb1(j,2)+sin(sum(Q1(j,:)))],Z2,'g'); hold on    
            %debug: this shows that length of the forearms vary according
            %to the time !!!!!! :(
%             disp(['length fore1 = ',num2str(norm(Elb1(j,:)-P1(j,:)))])
        delete(fore2) ; fore2 = patch('Vertices',V2{j,1},'Faces',templateForeF,'FaceColor','b','FaceAlpha',alpha); hold on
%         fore2 = plot3([Elb2(j,1);Elb2(j,1)+cos(sum(Q2(j,:)))],[Elb2(j,2);Elb2(j,2)+sin(sum(Q2(j,:)))],Z2,'b'); hold on      
            %debug:
%             disp(['length fore2 = ',num2str(norm(Elb2(j,:)-P2(j,:)))])
        delete(finger1) ; finger1 = patch('Vertices',V3{j,1},'Faces',templateForeF,'FaceColor','g','FaceAlpha',alpha); hold on
        delete(finger2) ; finger2 = patch('Vertices',V4{j,1},'Faces',templateForeF,'FaceColor','b','FaceAlpha',alpha); hold on
        
        delete(angle1_p1);
        angle1_p1 = plot_arc(0,Q1(j,1),orig1,r,c2_1/2,'g',scal);
        delete(angle2_p1);
        angle2_p1 = plot_arc(Q1(j,1),sum(Q1(j,:)),Elb1(j,:),r,c2_1/2,'g',scal);
        delete(angle1_p2);
        angle1_p2 = plot_arc(0,Q2(j,1),orig2,r,c2_2/2,'b',scal);
        delete(angle2_p2);
        angle2_p2 = plot_arc(Q2(j,1),sum(Q2(j,:)),Elb2(j,:),r,c2_2/2,'b',scal);
        
        delete(t1) ; delete(t2)
        t1 = text(ones(3,1)*(-0.25),[-0.00 -0.20 -0.40]',zeros(3,1), ...
                 {['q_{1_{player_{1}}}^{R_{0}}(t) = ',num2str(Q1(j,1)*180/pi),'°'], ...
                  ['q_{2_{player_{1}}}^{R_{0}}(t) = ',num2str(Q1(j,2)*180/pi),'°'], ...
                  ['|v|_{EE_{player_{1}}}^{R_{0}}(t) = ',num2str(normv1(j)),' m/s'] ...
                 }, ...
                'FontSize',txtftsize,'Color',[0 0.3 0]);
        t2 = text(ones(3,1)*0.6,[0.70 0.50 0.30]',zeros(3,1), ...
                 {['q_{1_{player_{2}}}^{R_{0}}(t) = ',num2str(Q2(j,1)*180/pi),'°'], ...
                  ['q_{2_{player_{2}}}^{R_{0}}(t) = ',num2str(Q2(j,2)*180/pi),'°'], ...
                  ['|v|_{EE_{player_{2}}}^{R_{0}}(t) = ',num2str(normv2(j)),' m/s'] ...
                 }, ...
                'FontSize',txtftsize,'Color','b');
        
        drawnow
%         pause(1)
%         if j==1 ; pause(20) ; end %debug
            if (j==i) ... %we arrive at the bumping moment
            && (collisionOrNo == 10 || occlusionOrNo == 21 || occlusionOrNo == 22)
                if EE1fore2;set(circle1,'Color','r'); set(fore2,'FaceColor','r'); end
                if EE2fore1;set(circle2,'Color','r'); set(fore1,'FaceColor','r'); end
                if EE1arm2;set(circle1,'Color','r'); set(arm2,'Color','r'); end
                if EE2arm1;set(circle2,'Color','r'); set(arm1,'Color','r'); end
                if foreArms;set(fore1,'FaceColor','r'); set(fore2,'FaceColor','r'); end
                if EEs; set(circle1,'Color','r'); set(circle2,'Color','r'); end
                if fore1arm2;set(fore1,'FaceColor','r'); set(arm2,'Color','r'); end
                if fore2arm1;set(fore2,'FaceColor','r'); set(arm1,'Color','r'); end
                if arms;set(arm1,'Color','r'); set(arm2,'Color','r'); end
                if fingers;set(finger1,'FaceColor','r'); set(finger2,'FaceColor','r'); end
                if fing1arm2;set(finger1,'FaceColor','r'); set(arm2,'Color','r'); end
                if fing2arm1;set(finger2,'FaceColor','r'); set(arm1,'Color','r'); end
                if fing1fore2;set(finger1,'FaceColor','r'); set(fore2,'FaceColor','r'); end
                if fing2fore1;set(finger2,'FaceColor','r');set(fore1,'FaceColor','r');end
                if vision1fing2;set(finger2,'FaceColor','r');set(vision1,'Color','r');end
                if vision2fing1;set(finger1,'FaceColor','r');set(vision2,'Color','r');end
                if vision1fore2;set(fore2,'FaceColor','r');set(vision1,'Color','r');end
                if vision2fore1;set(fore1,'FaceColor','r');set(vision2,'Color','r');end
                
                pause(1)
            end
    end

    plotTime = toc;
    disp(['The plot took ',num2str(plotTime-obs),' s.'])
end

yesorno = collisionOrNo + occlusionOrNo;
end