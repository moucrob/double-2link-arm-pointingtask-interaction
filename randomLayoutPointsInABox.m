function [posMatNx2] = randomLayoutPointsInABox(N, sizeBoxX, sizeBoxY)
%For generate in a uniformly RANDOM ARRANGED WAY N points in a box :

%stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly/50746409#50746409
clc
nbPoints = 1000;

set1 = sizeBoxX*rand(nbPoints,1);
set2 = sizeBoxY*rand(nbPoints,1); %each time we recall rand it generates new set

% ampl = 1;
% 
% thetas = 2*pi*set2;
% 
% %For a circle :
% R = 1;

%% DOESN'T WORK
% %For a square :
% R = [];
% for i=1:numel(thetas)
% %     disp(['i = ',num2str(i)])
%     between0and8 = thetas(i)/(pi/4);
%     borneinf = floor(between0and8);
% %     bornesup = ceil(between0and8);
%     if any([0 3 4 7] == borneinf)
%         R(end+1,1) = sqrt( 1 + (sin(thetas(i)))^2 );
%     else
%         R(end+1,1) = sqrt( 1 + (cos(thetas(i)))^2 );
%     end
% end
% r = ampl*R.*sqrt(set1);
%%
%%
% x = r.*cos(thetas);
% y = r.*sin(thetas);
% figure
% for i=1:numel(x)
%     plot(x,y,'.') ; hold on
% end

% plot(set1,set2,'.')
posMatNx2 = [set1 set2];


% th=0:(2*pi/1000):2*pi;
% xline = R*cos(th);
% yline = R*sin(th);
% plot(xline,yline)