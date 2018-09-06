set(0,'DefaultFigureWindowStyle','docked')
figure
set(0,'DefaultFigureWindowStyle','normal')

% test = [[1;0;0;1]];
test = [[2;0;0;1]];
p = T(pi/2,[3;3;0])*test;

plt = plot3(test(1,:),test(2,:),test(3,:),'k.') ; hold on
pl = plot3(p(1,:),p(2,:),p(3,:),'r.');
xlim([0 5]) ; ylim([0 5]) ; zlim([0 inf]) ; grid on ; grid minor
hold on

%So in summary, T(angle,pos) translate the initial attached (to 0) vector
%in the position pos (w.r.t the global frame 0), and then rotate it!

%So if we want t to be a "sliding translation" vector, then we have to add
%in input of T the current pos of the vector to be translated!