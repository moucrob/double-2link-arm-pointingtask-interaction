%{
Database
    %Rectangle polyhedron :
    %V = [1 1 0 ; 1 3 0 ; 0 3 0 ; 0 1 0 ; 1 1 1 ; 1 3 1 ; 0 3 1 ; 0 1 1];
    %F = [1 2 6 5 ; 2 3 7 6 ; 3 4 8 7 ; 4 1 5 8 ; 1 2 3 4 ; 5 6 7 8];

    %Octahedron :
    %V = [-1 0 0 ; 0 1 0 ; 1 0 0 ; 0 -1 0 ; 0 0 -1 ; 0 0 1];
    %F = [ 1 2 5 ; 1 2 6 ; 2 3 5 ; 2 3 6 ; 3 4 5 ; 3 4 6 ; 1 4 5 ; 1 4 6];
%}
%About the ensuring of the outwarding normal, it implies that POLYHEDRON
%MUST BE CONVEX (see my A4 drawing normally included in my notebook)
function Nx3 = outwardingFaceNormalVectorsOfAConvexPolyhedron(faces,vertices,display)
%because of the use of patch() for my work, I assume all faces of the
%polyhedron have the same number of vertices.

clc

%For the computation of the normal:
if size(vertices,2) ~= 3
    disp('ERROR : Cross-product in Matlab must have in input 3D arrays')
    return
end

%(For the ensuring of good orientation of the normal:)
%Establishing of the list of all vertices:
list = 1:size(vertices,1);
    
%For validating the good behaviour of the function:
if display == 1
    nbVertices = size(faces,2);
    alpha = 0.5; %transparency
    set(0,'DefaultFigureWindowStyle','docked')
    figure ; grid on ; grid minor
    set(0,'DefaultFigureWindowStyle','normal')
    set(gca,'DataAspectRatio',[1 1 1]) ; view(3)
    xlabel('x') ; ylabel('y') ; zlabel('z')
    patch('vertices',vertices,'faces',faces,'FaceColor','g','FaceAlpha',alpha);
    hold on
    uselessplot = plot3(0,0,0,'.') ; tmp = uselessplot;
end

for i=1:size(faces,1) %for each face,
    %% Compute its normal:
    %(One face has at least 3 vertices)
    a1 = vertices( faces(i,2) , : ) - vertices( faces(i,1) , : );
    a2 = vertices( faces(i,3) , : ) - vertices( faces(i,2) , : );
    nOld = cross(a1,a2);
    %% Ensure it is the the outwarding and not inwarding one:
    % math.stackexchange.com/questions/183030/given-a-tetrahedron-how-to-find-the-outward-surface-normals-for-each-side
    %Goal : pick one vertice of the convex polyhedron not belonging to this
    %face.
    remainingSet = list ; remainingSet([faces(i,:)]) = [];
    v_endpoint = vertices(remainingSet(randperm(numel(remainingSet),1)),:);
    v = v_endpoint - vertices(faces(i,1),:);
%     disp(['sign(dot(nOld,v)) = ',num2str(sign(dot(nOld,v)))])
    if sign(dot(nOld,v)) < 0
        nNew = nOld;
    else %dot>0
        nNew = -nOld;
    end
    Nx3(i,:) = nNew;
    
    %For validating the good behaviour of the function:
    if display == 1
        midPoint = sum(vertices( faces(i,:) , : ),1)/nbVertices;
        quiver3(midPoint(1),midPoint(2),midPoint(3),nNew(1),nNew(2),nNew(3),'r-','LineWidth',2,'AutoScale','off')
        
        delete(tmp) ; tmp = [];
        tmp(end+1) = quiver3(vertices(faces(i,1),1),vertices(faces(i,1),2),vertices(faces(i,1),3),a1(1),a1(2),a1(3),'b-','LineWidth',2,'AutoScale','off');
        tmp(end+1) = quiver3(vertices(faces(i,2),1),vertices(faces(i,2),2),vertices(faces(i,2),3),a2(1),a2(2),a2(3),'b-','LineWidth',2,'AutoScale','off');
        tmp(end+1) = quiver3(vertices(faces(i,1),1),vertices(faces(i,1),2),vertices(faces(i,1),3),nOld(1),nOld(2),nOld(3),'r:','LineWidth',2,'AutoScale','off');
        tmp(end+1) = quiver3(vertices(faces(i,1),1),vertices(faces(i,1),2),vertices(faces(i,1),3),v(1),v(2),v(3),'k:','LineWidth',2,'AutoScale','off');
        
        drawnow ; hold on ; %pause()
    end
end