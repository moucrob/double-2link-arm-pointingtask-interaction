%useless
figure ; set(gca,'DataAspectRatio',[1 1 1]) ; grid on ; grid minor
vFo = [b2_1 0;
       b2_1 lenFore;
       -b2_1 lenFore;
       -b2_1 0];
face = [1 2 3 4];
vFi = [-b2_1 lenFore] + ...
      [b3_1 0;
      b3_1 projectedLenFing;
      -b3_1 projectedLenFing;
      -b3_1 0];
patch('vertices',vFo,'faces',face,'FaceColor','c') ; hold on
patch('vertices',vFi,'faces',face,'FaceColor','c') ; hold on
plot([0 0],[0 lenFore],'k-') ; hold on
hypo_x12 = [0 -b2_1] ; hypo_y12 = [0 lenFore+projectedLenFing];
plot(hypo_x12,hypo_y12,'k-') ; hold on
[epsilon,see] = perpendicularDecenteringEE(b2_1,lenFore+projectedLenFing);
disp(['Length of the shortest distance between O_2 and the EE decentered = ',num2str(see)])
hypo_conc = [hypo_x12;
             hypo_y12];
hypo_vec = hypo_conc(:,2) - hypo_conc(:,1) ; normhypo = norm(hypo_vec);
disp(['normhypo = ',num2str(normhypo)])
plot_arc(pi/2,(pi/2)+epsilon,[0 0],4*lenFore/5,0,'r',lenFore/20)