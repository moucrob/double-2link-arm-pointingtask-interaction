function [angle,a1_real,projectedLenFing] = whichDesignAngleRotY(c2,a3,c3,distFromElbToEE,demo)
%See my notebook p. 12 July for understand with the drawing
%a3 = length of the finger
%c3 = height of the finger
%c2 = height forearm (more precisely of the block hand)
hypo = sqrt(a3^2 + (c3/2)^2); %(hypothenuse)

epsilon = asin((c3/2)/hypo);

%The first calculus line is:
%hypo*sin(pi+eps+x?)= -c2;
angle = asin(c2/hypo)-epsilon;

%Hence:
projectedLenFing = hypo*cos(angle+epsilon);
a1_real = distFromElbToEE - projectedLenFing;

%% DEMO
if demo
    figure ; alpha = 0.5 ; set(gca,'DataAspectRatio',[1,1,1])
    ForeV = [ 0   0;
              0  -c2;
             0.4 -c2; %0.4 is not accurate
             0.4  0];
    ForeF = [1 2 3 4];
    patch('Vertices',ForeV,'Faces',ForeF,'FaceColor','c','FaceAlpha',alpha) ; hold on

    HandV = [ 0   0;
              0  -c2;
             0.1 -c2; %0.4 is not accurate
             0.1  0];
    patch('Vertices',HandV,'Faces',ForeF,'FaceColor','c','FaceAlpha',alpha) ; hold on

    FingV = [ 0   c3/2;
             -a3  c3/2;
             -a3 -c3/2;
              0  -c3/2];
    patch('Vertices',FingV,'Faces',ForeF,'FaceColor','c','FaceAlpha',alpha) ; hold on
    toSee = 2 ; plot(toSee*[0;-a3],toSee*[0;-c3/2],'k--') ; hold on %hypo
    plot(toSee*[0;-a3],toSee*[0;0],'k--') ; hold on %bone
    plot_arc(pi,pi+epsilon,[0 0],(7/4)*a3,0,'b',a3/5)

    plot_arc(pi+epsilon,pi+epsilon+angle,[0 0],hypo,0,'r',hypo/5)
    grid on ; grid minor ; xlabel('x') ; ylabel('y') ; zlabel('z')

    disp(['angle = ',num2str(angle*180/pi),'°'])
end
end