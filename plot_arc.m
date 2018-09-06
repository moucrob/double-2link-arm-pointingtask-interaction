function P = plot_arc(qini,qend,center,r,alti,color,scal)
%centered on [cX;cY], of radius r, starting from qini to qend
% actually useless afterward:
% https://fr.mathworks.com/help/matlab/ref/matlab.graphics.shape.arrow-properties.html
cX = center(1) ; cY = center(2);

s = linspace(qini,qend,150)';
% global x y %just for testing something, to remove!
x = cX + r*cos(s);
y = cY + r*sin(s);
z = alti*ones(numel(s),1);

R = @(th) ...
    [cos(th) -sin(th) 0;
     sin(th)  cos(th) 0;
        0        0    1];
slipperyPreArrowHead = [x(end-1)-x(end) ; y(end-1)-y(end) ; 0]; %z = 0 before scaling, otherwise the scale rises the segment!!
%one portion among the 150 portions of the arc might be more or less long
%depending of the space between qini and qend, hence arrowheads depend on
%delta_q, thing that I don't want since I brought the parameter scal, so I
%have to normalize this vec:
slipperyPreArrowHead = slipperyPreArrowHead/norm(slipperyPreArrowHead);
slipperyPreArrowHead_above = R(pi/4)*scal*slipperyPreArrowHead;
slipperyPreArrowHead_under = R(-pi/4)*scal*slipperyPreArrowHead;
ending = [x(end) y(end) alti]';
preArrowHead_above = slipperyPreArrowHead_above + ending;
preArrowHead_under = slipperyPreArrowHead_under + ending;
arrowHead = [preArrowHead_above' ; ...
             ending' ; ...
             preArrowHead_under'];
         

P = plot3(x,y,z,color,arrowHead(:,1),arrowHead(:,2),arrowHead(:,3),color); hold on
% plot3()% ; hold on

%[x_ini x_end] for orienting the arrowhead :
% X_DS = [x(end-1) x(end)]; %Where we want the annotation ON THE PLOT
% Y_DS = [y(end-1) y(end)]; %(DS stands for data space)

%annotation function requires positions but in a strange scaled way :
%bottom left of the WINDOW (~= plot) is considered as (0;0) while top right
%is considered as (1;1):
% [X_NFU,Y_NFU] = ds2nfu(X_DS,Y_DS); %(NFU stands for normalized figure units)

% annotation('arrow', 'Color', color, ...
%            'LineStyle','none', ...
%            'HeadStyle','cback3', ...
%            'HeadLength',head_size, 'HeadWidth',head_size, ...
%            'Units','normalized', ...
%            'X',X_NFU,'Y',Y_NFU)


% position = [x(end-1) y(end-1) x(end) y(end)];
% h = annotation('arrow'); 
% set(h,'parent', gca, 'Position', position, ...
%     'HeadLength',head_size, 'HeadWidth', head_size, ...
%     'HeadStyle',head_style, 'LineStyle','none', 'Color',color);

% quiver3(x(end-1), y(end-1), 0, ... %starting point
%         x(end)-x(end-1), y(end)-y(end-1), 0, ... %amplitudes
%         'Color',color, 'AutoScale','off', 'MaxHeadSize',head_size) %otherwise it's miss-scaled


