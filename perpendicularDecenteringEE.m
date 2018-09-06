function [DqLast,lenTmpSkeleton] = perpendicularDecenteringEE(shiftDist,lenCentralSkeleton)
%{
When you want to know from how many radians you should adjust your verry
%last link if you want its end effector makes a .--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    |
%instead of just a straight line:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    .
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    |
%Then one should compute the IK with as last link the hypothenuse:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% .--.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  \ |
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   \|
%with then rotating the central skeleton by q2 wrongly computed
%more or less Delta q_last (here minus if decentering to the left and if we
%use normally oriented angles starting from the trigonometric circle frame
%[x=1,y=0].
%}

hypo = sqrt(shiftDist^2 + lenCentralSkeleton^2);

%utilisation of Al-Kashi (Generalized Pythagore) Theorem knowing 3edges in
%order to determine the angle:
DqLast = acos( (hypo^2 + lenCentralSkeleton^2 - shiftDist^2)/(2*hypo*lenCentralSkeleton) );
lenTmpSkeleton = hypo;
end