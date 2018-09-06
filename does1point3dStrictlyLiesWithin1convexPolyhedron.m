%{
To try : initialize :
%Octahedron :
    %V = [-1 0 0 ; 0 1 0 ; 1 0 0 ; 0 -1 0 ; 0 0 -1 ; 0 0 1];
    %F = [ 1 2 5 ; 1 2 6 ; 2 3 5 ; 2 3 6 ; 3 4 5 ; 3 4 6 ; 1 4 5 ; 1 4 6];
and ask:
    yesorno = does1point3dStrictlyLiesWithin1convexPolyhedron([-1.5 -1.5 0],F,V,1)
%}
function yesorno = does1point3dStrictlyLiesWithin1convexPolyhedron(pointOfInterest,faces,vertices,disp)
disp = disp; %displays whether a point lies within
display = 0; %displays outwarding face normal vectors

%In order to validate the right behaviour of the function:
if disp == 1
    nbVertices = size(faces,2);
    alpha = 0.5; %transparency
    set(0,'DefaultFigureWindowStyle','docked')
    figure ; grid on ; grid minor
    set(0,'DefaultFigureWindowStyle','normal')
    set(gca,'DataAspectRatio',[1 1 1]) ; view(3) ; camorbit(20,0,'data',[0 0 1])
    xlabel('x') ; ylabel('y') ; zlabel('z')
    patch('vertices',vertices,'faces',faces,'FaceColor','g','FaceAlpha',alpha);
    hold on
    uselessplot = plot3(0,0,0,'.') ; tmp = uselessplot;
end

Nx3 = outwardingFaceNormalVectorsOfAConvexPolyhedron(faces,vertices,display);
%(Each normal (row of Nx3) index is the same as the one of the Face (row of
%faces).
if size(pointOfInterest,1) ~= 1 ; transpose(pointOfInterest) ; end
yesorno = 1; %First let's assume the point strict-lies inside.
for i=1:size(faces,1)
    if disp == 1
        midPoint = sum(vertices( faces(i,:) , : ),1)/nbVertices; %optionnal
        %line, but it's less confusing when visualizing
        faceToPointOfInterest = pointOfInterest - midPoint;
        delete(tmp) ; tmp=[];
        tmp(end+1) = quiver3(midPoint(1),midPoint(2),midPoint(3),faceToPointOfInterest(1),faceToPointOfInterest(2),faceToPointOfInterest(3),'k:','LineWidth',2,'AutoScale','off');
        tmp(end+1) = quiver3(midPoint(1),midPoint(2),midPoint(3),Nx3(i,1),Nx3(i,2),Nx3(i,3),'r-','LineWidth',2,'AutoScale','off');
    else
        faceToPointOfInterest = pointOfInterest - vertices( faces(i,1) , : );
    end
    if sign(dot(Nx3(i,:),faceToPointOfInterest)) > 0 
        yesorno = 0;
        return
    end
    if disp == 1 ; pause() ; end
end

end