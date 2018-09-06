function Q = IKright2RarmWithoutAlKashi(EEpos,upperArmLen,foreArmLen)
%www.hesmer.org/uploads/RobotArm/Inverse%2520Kinematics%2520for%2520Robot%2520Arm.pdf
%instead of this method I was told in course :
%https://robotacademy.net.au/lesson/inverse-kinematics-for-a-2-joint-robot-arm-using-geometry/

x = EEpos(1) ; y = EEpos(2) ; a1 = upperArmLen ; a2 = foreArmLen;
% disp(['x = ', num2str(x)])
% disp(['y = ', num2str(y)])

interm = (x^2 + y^2 - a1^2 - a2^2)/(2*a1*a2);
%when (x;y) is located within a circle centered on the origin and with a
%small radius, then a1² > x² and a2² > y² so interm remains negative.
%Now let's see if interm² > 1...
% disp(['interm = ', num2str(interm)])
% disp(['r = ', num2str(sqrt(x^2 + y^2))])
% disp(['a1 + a2 = ', num2str(a1+a2)])

q2 = atan2(sqrt(1-interm^2) , interm); % +sqrt and not -sqrt bcuz right arm and not left.
q1 = atan2(y,x) - atan2(a2*sin(q2) , a1+a2*cos(q2));
Q = [q1 q2];

end