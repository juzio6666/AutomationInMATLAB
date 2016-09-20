plant = servo();

umin = -220;
umax = 220;
init_state = [0, 0, 0, 0];

% simulation parameters
stages = 20;
stage_time = 100;
stage_resolution = 10;
T = 0:1/stage_resolution:(stage_time*stages-1/stage_resolution);
y = [];
t = 0:1/stage_resolution:(stage_time-1/stage_resolution); 
x = init_state;
Y = [];
U = [];

% simulation
for i=1:stages
    r = umin + (umax-umin).*rand(1,1);
    [y,~,x]=lsim(plant,r*ones(size(t)),t,x(size(x,1),:));
    Y = [Y; y];
    U = [U, r*ones(size(t))];
end

% results
figure;
subplot(2,1,1);
plot(T,Y)
xlim([0, max(T)]);
subplot(2,1,2);
plot(T,U)
xlim([0, max(T)]);
ylim([umin, umax]);

% step response
[y,~,x]=lsim(plant,1*ones(size(t)),t,[0,0,0,0]);
figure;
plot(t,y)