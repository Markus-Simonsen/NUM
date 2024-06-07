[X,Y] = meshgrid(-1:0.1:1,-1:0.1:1);
Z = X.^4 - 2 * Y.^4;

v = readvars('V.csv');
u = readvars('U.csv');
r = v.^4 - 2 * u.^4; 
view(3)
hold on;
s = surf(X, Y, Z,'FaceAlpha',0.5);
plot3(v,u,r,'Color','b','linewidth',3);
xlabel('v');
ylabel('u');
zlabel('r');
