clc

alpha = 0.5; %transparency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TO DO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterationsAllowed = 2; %IMPORTANT PARAMETER TO STUDY 

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
        
x2 = [3 ; 4]; %-1; %-2.5;
%-1 implies that one vertex of the segment can exactly touch a face and be
%detected!
y2 = [0.5 ; 0.5];
z2 = [0.5 ; 0.5];
V2 = [x2,y2,z2];
F2 = [1 2];
        
set(0,'DefaultFigureWindowStyle','docked')
figure ; grid on ; grid minor
set(0,'DefaultFigureWindowStyle','normal')
set(gca,'DataAspectRatio',[1 1 1]) ; view(3)
xlabel('x') ; ylabel('y') ; zlabel('z')

S1obj = patch('Vertices',V1,'Faces',F1,'FaceColor','g','FaceAlpha',alpha);
hold on
S2obj = patch('Vertices',V2,'Faces',F2,'FaceColor','k');
hold on

collisionFlag12 = GJK(S1obj,S2obj,iterationsAllowed)

x3 = repmat([0.5 ; 1.5 ; 1.5 ; 0.5],2,1); %+1; %repeat the vector. Putting ,2,1 or ,1,2
%is quite unpredictable, didnt understood the logic, so I prefer not comment this way...
y3 = repmat([0.3 ; 0.3 ; 0.7 ; 0.7],2,1);
z3 = [0.3*ones(4,1);0.7*ones(4,1)];
V3 = [x3,y3,z3];
F3 = [1 2 6 5 ; ... %Face1 connects vertices (rows of V1) number 1,2,5,6
            2 3 7 6 ; ... 
            3 4 8 7 ; ...
            4 1 5 8 ; ...
            1 2 3 4 ; ...
            5 6 7 8];
S3obj = patch('Vertices',V3,'Faces',F3,'FaceColor','b','FaceAlpha',alpha);
collisionFlag13 = GJK(S1obj,S3obj,iterationsAllowed)

%Conclusion: case dealing with segment lying within the polyhedron still
%has to be detected!