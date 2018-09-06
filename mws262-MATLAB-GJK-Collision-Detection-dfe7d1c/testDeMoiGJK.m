clc
% https://fr.mathworks.com/help/matlab/visualize/multifaced-patches.html

alpha = 0.5; %transparency

x1 = repmat([0 ; 2 ; 2 ; 0],2,1); %repeat the vector. Putting ,2,1 or ,1,2
%is quite unpredictable, didnt understood the logic, so I prefer not comment this way...
y1 = repmat([0 ; 0 ; 1 ; 1],2,1);
z1 = [zeros(4,1);ones(4,1)];
V1 = [x1,y1,z1];
F1 = [1 2 6 5 ; ... %Face1 connects vertices (rows of V1) number 1,2,5,6
            2 3 7 6 ; ... 
            3 4 8 7 ; ...
            4 1 5 8 ; ...
            1 2 3 4 ; ...
            5 6 7 8];
origBlock1_0 = V1(1,:)'

% x2 = x1+4;
% x2 = x1+2;
% y2 = y1; ; z2 = z1;
% V2 = [x2,y2,z2];

t3x1 = [1;1;0];
T = @(theta,origInFrame0,t3x1) ...
    [cos(theta) -sin(theta) 0 t3x1(1); %origInFrame0(1)+t3x1(1);
     sin(theta) cos(theta)  0 t3x1(2); %origInFrame0(2)+t3x1(2);
          0          0      1 t3x1(3); %origInFrame0(3)+t3x1(3);
          0          0      0    1   ];
toTrim = T(pi/2,origBlock1_0,t3x1)*[V1';ones(1,size(V1,1))];
% V2 = toTrim(1:3,:)';
V1 = toTrim(1:3,:)';
                   
F2 = F1;

set(0,'DefaultFigureWindowStyle','docked')
figure ; grid on ; grid minor
set(0,'DefaultFigureWindowStyle','normal')
set(gca,'DataAspectRatio',[1 1 1]) ; view(3)
xlabel('x') ; ylabel('y') ; zlabel('z')
%Only for this snippet:
% xlim([0 inf]) ; ylim([0 inf])

S1obj = patch('Vertices',V1,'Faces',F1,'FaceColor','g','FaceAlpha',alpha);
% S2obj = patch('Vertices',V2,'Faces',F2,'FaceColor','b','FaceAlpha',alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TO DO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterationsAllowed = 2; %IMPORTANT PARAMETER TO STUDY 

% tic
% collisionFlag = GJK(S1obj,S2obj,iterationsAllowed);
% t = toc ; disp(['The collision function took ', num2str(t),' sec.'])
% disp(['collisionFlag = ', num2str(collisionFlag)])