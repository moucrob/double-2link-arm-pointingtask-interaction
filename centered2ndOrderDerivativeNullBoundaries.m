function der = centered2ndOrderDerivativeNullBoundaries(pos,dt)
%Centered derivatives of order 2 (more accurate than backward or forward of
%order 1). With added null boundaries for the realness.
%(dt is a scalar fixed)

sizepos = numel(pos);

der = [0];
if sizepos > 2
    for i=2:sizepos-1 ; der(i,1) = (pos(i+1)-pos(i-1))/2*dt ; end
else
    disp('Be careful, pos is not at least 3-frames-long, hence it is not possible to use central derivatives')
end
der(end+1) = 0;

end