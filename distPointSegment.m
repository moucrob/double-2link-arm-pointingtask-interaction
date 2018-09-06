%% Credits to:
%bretus :
%https://developpez.net/forums/d1304134/general-developpement/algorithme-mathematiques/mathematiques/distance-point-segment/
%% Function:
function len = distPointSegment(argx4,point)
xA = argx4(1) ; yA = argx4(2); A = [xA,yA];
xB = argx4(3) ; yB = argx4(4); B = [xB,yB];
xP = point(1) ; yP = point(2); P = [xP,yP];

AP = P-A ; AB = B-A;
%     %debug
%     p(1) = quiver(xA,yA,AB(1),AB(2),'AutoScale','off') ; hold on
%     p(2) = quiver(xA,yA,AP(1),AP(2),'AutoScale','off') ; hold on

%TO NORMALIZE AB IS, AFTERWARDS, IMPORTANT ! ...
s = dot(AP,AB/norm(AB)); %distance pointing from A to the orthogonal projection of P, P', on AB
%(if s < 0, P' doesn't belong to AB)
%(if 0<s<|AB|, P' does.)
%     %debug
%     disp(['s = ',num2str(s)])
    

s = s/norm(AB); %(if 0<s<1, P' does.)
%     %debug
%     disp(['s = ',num2str(s)])

s = min( max(0,s) , 1 );
%     %debug
%     disp(['s = ',num2str(s)])

Pp = A + AB*s; %Pp stands for P'
%     %debug
%     p(3) = plot(Pp(1),Pp(2),'.','MarkerSize',20) ; hold on
%     PPp = Pp-P;
%     p(4) = quiver(xP,yP,PPp(1),PPp(2),'k:','AutoScale','off') ; hold on

len = norm(P-Pp);
%     %debug
%     disp(["|PP'| = ",num2str(len)])
end