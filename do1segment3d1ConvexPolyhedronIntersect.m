function yesorno = do1segment3d1ConvexPolyhedronIntersect(segCoord2x3,vertices,faces,plt)

F1 = [1 2];

S1obj = patch('Vertices',segCoord2x3,'Faces',F1,'visible','off');
S2obj = patch('Vertices',vertices,'Faces',faces,'visible','off'); %,'FaceColor','g','FaceAlpha',alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   TO DO    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iterationsAllowed = 2; %IMPORTANT PARAMETER TO STUDY 

%Does the segment penetrates and is stucked into one of the faces?
yesorno = GJK(S1obj,S2obj,iterationsAllowed);

if yesorno == 0
    %Does the segments lies within the convex Polyhedron?
    %i.e do one of the 2 extremity points lie within the CP?
    disp = 0;
    yesorno = does1point3dStrictlyLiesWithin1convexPolyhedron(segCoord2x3(1,:),faces,vertices,disp);
end