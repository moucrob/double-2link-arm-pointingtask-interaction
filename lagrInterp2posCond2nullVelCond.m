function P = lagrInterp2posCond2nullVelCond(posx2,duration,t)
%interpolate a trajectory between 2 boundary positions given 2 null
%velocity conditions.
dim = numel(posx2)/2;
posIni = posx2(1:dim); posFin = posx2(dim+1:end);
% disp(['posIni = ',num2str(posIni)]) ; disp(['posFin = ',num2str(posFin)])

if t<=duration
%     function I found for stephane caron assignment (big big mistake on the
%     last row of my matrix T.a = cond :'(
%     P =   (1/duration^3)*(posFin - (2/3)*posIni)*t^3 ...
%         - (1/(3*duration^2))*posIni*t^2 + posIni;
    P = (-(2/duration^3)*t^3+(3/duration^2)*t^2)*(posFin-posIni) + posIni;
%     disp(['P = ',num2str(P)])
else
    disp('Error, t > duration')
end

end