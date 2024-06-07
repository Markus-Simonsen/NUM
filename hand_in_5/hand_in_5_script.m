

u = -1:0.1:1;
v = [-1:0.1:1]';
r_func = @(u,v) v.^4-2.*u.^4;
r = r_func(u,v)
disp(r)

f1 = figure("NumberTitle","off");
f1.Name = "Figure in 3D"
figure(f1);
surf(u,v,r)
xlabel("u")
ylabel("v")
zlabel("r")
hold on
plot3(-0.9,-0.85, r_func(-0.9, -0.85), 'or','LineWidth',3)
plot3(0.8,-0.9, r_func(0.8, -0.9), 'or','LineWidth',3)
hold off


%%
a = -0.9;
alpha = -0.85;
b = 0.8;
beta = -0.9

umesh = linspace(-0.9, 0.8, 10);
solinit = bvpinit(umesh, [0; 1]);

sol4c = bvp5c(@bvpfcn, @bcfcn, solinit)

f2 = figure;
figure(f2)
plot(sol4c.x,sol4c.y(1,:))
%legend('1','2')

figure(f1);
hold on
plot3(sol4c.x,sol4c.y(1,:),r_func(sol4c.x,sol4c.y(1,:)), '*')
hold off
%%
function res = bcfcn(Va,Vb)
res = [Va(1)+0.85
       Vb(1)+0.9];
end
%%
function dVdu = bvpfcn(u,V)
dVdu = [V(2)
       48*(V(1)^3+2*u^3*V(2))*(2*u^2-V(1)^2*V(2)^2)];
end