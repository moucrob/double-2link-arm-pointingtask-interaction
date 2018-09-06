% Ry = @(b) [cos(b) 0 sin(b) ; 0 1 0 ; -sin(b) 0 cos(b)];
% Rz = @(a) [cos(a) -sin(a) 0 ; sin(a) cos(a) 0 ; 0 0 1];
% 
% extr1 = [0;0;0]; extr2 = [1;0;0]; %bones of the template initial polyhedrons
% toPlot = [extr1';extr2'];
% figure ; plot3(toPlot(:,1),toPlot(:,2),toPlot(:,3),'r') ; hold on
grid on ; grid minor ; xlabel('x'), ylabel('y') ; zlabel('z')
% 
% first = pi/4;
% extr1 = Ry(first)*extr1 ; extr2 = Ry(first)*extr2 ; hold on
% toPlot = [extr1';extr2'];
% plot3(toPlot(:,1),toPlot(:,2),toPlot(:,3),'g')
% 
% second = pi/2;
% extr1 = Rz(second)*extr1 ; extr2 = Rz(second)*extr2;
% butBefore = [3;3;1];
% extr1 = extr1 + butBefore ; extr2 = extr2 + butBefore;
% toPlot = [extr1';extr2'];
% plot3(toPlot(:,1),toPlot(:,2),toPlot(:,3),'b') ; hold on

rx2 = [1 1]*0.005;
dimsArms = [ones(2,1)*0.32 ...
            ones(2,1)*0.4 ones(2,1)*0.075 ones(2,1)*0.06 ...
            ones(2,1)*0.1 ones(2,2)*2*rx2(1)*cos(pi/4) ]; %(m) (logic below)
a3_1 = dimsArms(1,5); %length 1st finger
b3_1 = dimsArms(1,6); %width 
c3_1 = dimsArms(1,7); %height

templateFingX = repmat(a3_1*[0 ; 1 ; 1 ; 0],2,1); %repeat the vector. Putting ,2,1 or ,1,2
%is quite unpredictable, didnt understood the logic, so I prefer not comment this way...
templateFingY = repmat((b3_1/2)*[-1 ; -1 ; 1 ; 1],2,1);
templateFingZ = c3_1*[zeros(4,1);ones(4,1)];
templateFingV = [templateFingX,templateFingY,templateFingZ];
templateFingF = [1 2 6 5 ; ... %Face1 connects vertices (rows of V1) number 1,2,5,6
                 2 3 7 6 ; ... 
                 3 4 8 7 ; ...
                 4 1 5 8 ; ...
                 1 2 3 4 ; ...
                 5 6 7 8 ];
patch('Vertices',templateFingV,'Faces',templateFingF,'FaceColor','m') ; hold on

Ry = @(phi) ...
    [cos(phi) 0 sin(phi); ...
        0     1    0    ; ...
    -sin(phi) 0 cos(phi)];
whichAngle = pi/4;
templateFingV = Ry(whichAngle)*templateFingV' ; templateFingV = templateFingV';
patch('Vertices',templateFingV,'Faces',templateFingF,'FaceColor','m') ; hold on

toBeTransformed = [templateFingV';ones(1,size(templateFingV,1))];
Tz = @(theta,t3x1) ...
    [cos(theta) -sin(theta) 0 t3x1(1);
     sin(theta) cos(theta)  0 t3x1(2);
          0          0      1 t3x1(3);
          0          0      0    1   ];
theta = pi/2 ; t3x1 = [0.3;0.3;0.1];
templateFingV = Tz(theta,t3x1)*toBeTransformed;
templateFingV(4,:) = [] ; templateFingV = templateFingV';
patch('Vertices',templateFingV,'Faces',templateFingF,'FaceColor','m')
