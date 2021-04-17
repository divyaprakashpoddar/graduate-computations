clear; clc; close;

x = @(t) sin(3*t);
y = @(t) sin(2*t);
z = @(t) sin(2*t + pi/2);

t = 0:0.1:600;

X = x(t); Y = y(t); Z = z(t);


Rx = @(a) [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
Ry = @(b) [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
Rz = @(g) [cos(g) -sin(g) 0; sin(g) cos(g) 0; 0 0 1];

Po = [X' Y' Z'];
Po = Po*Rz(pi/2);

figure(1)
plot3(Po(:,1),Po(:,2),Po(:,3),'.')
xlabel('x'); ylabel('y'); zlabel('z')

for i = -10*pi:0.05:10*pi
    P = Po*Ry(i);
    plot3(P(:,1),P(:,2),P(:,3),'.')
    view(2)
    axis off
    drawnow
end



