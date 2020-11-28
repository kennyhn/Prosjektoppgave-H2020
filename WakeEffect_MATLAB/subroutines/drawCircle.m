function [circleX, circleY, pathX, pathY] = drawCircle(posX,posY,scalex,scaley)

thetaMin = -pi;
thetaMax = pi;
r = 1;
r_path = 1.1;
dtheta = 0.01;

theta_vec = thetaMin:dtheta:thetaMax;
r_vec = r*ones(1,length(theta_vec));
pathIndx = 1:length(r_vec);
r_path_vec = r_path*ones(1,length(pathIndx));

[circleX, circleY] = pol2cart(theta_vec, r_vec);
circleX = circleX*scalex + posX;
circleY = circleY*scaley + posY;
[pathX, pathY] = pol2cart(theta_vec(pathIndx), r_path_vec);
pathX = pathX*scalex + posX;
pathY = pathY*scaley + posY;

end