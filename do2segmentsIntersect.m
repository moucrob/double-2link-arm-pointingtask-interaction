function yesorno = do2segmentsIntersect(argx8)

x1start = argx8(1);
y1start = argx8(2);
x2start = argx8(3);
y2start = argx8(4);
x1end = argx8(5);
y1end = argx8(6);
x2end = argx8(7);
y2end = argx8(8);
                                    
A = [[1 1 1]; [x1start x1end x2start]; [y1start y1end y2start]];

B = [1 1 1;
     x1start x1end x2end;
     y1start y1end y2end];

det1 = det(A)*det(B);

C = [1 1 1;
     x1start x2start x2end;
     y1start y2start y2end];

D = [1 1 1;
     x1end x2start x2end;
     y1end y2start y2end];

det2 = det(C)*det(D);

if (det1 <=0) && (det2 <=0)
    yesorno = 1;
else
    yesorno = 0;
end

end