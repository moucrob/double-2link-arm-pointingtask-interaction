function positions = circle2D(x,y,r,nbIntervals)

th = 0 : 2*pi/nbIntervals : 2*pi;
positions(:,1) = r*cos(th) + x;
positions(:,2) = r*sin(th) + y;

end