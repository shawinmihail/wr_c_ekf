% clc
% clear
close all

%% C
path = 'est_state_log.csv';
data = csvread(path);

x = data(:,1);
y = data(:,2);
z = data(:,3);

vx = data(:,4);
vy = data(:,5);
vz = data(:,6);

dx1 = data(:,7);
dy1 = data(:,8);
dz1 = data(:,9);

dx2 = data(:,10);
dy2 = data(:,11);
dz2 = data(:,12);

mx = data(:,13);
my = data(:,14);
mz = data(:,15);

mvx = data(:,16);
mvy = data(:,17);
mvz = data(:,18);

mdx1 = data(:,19);
mdy1 = data(:,20);
mdz1 = data(:,21);

mdx2 = data(:,22);
mdy2 = data(:,23);
mdz2 = data(:,24);


%% PLOT
figure
hold on
plot(mx, 'r')
% plot(my, 'g')
% plot(mz, 'b')
plot(x, 'r--')
% plot(y, 'g--')
% plot(z, 'b--')

figure
hold on
plot(mvx, 'r')
% plot(mvy, 'g')
% plot(mvz, 'b')
plot(vx, 'r--')
% plot(vy, 'g--')
% plot(vz, 'b--')


figure
hold on
plot(mdx1, 'r')
% plot(mdy1, 'g')
% plot(mdz1, 'b')
plot(dx1, 'r--')
% plot(dy1, 'g--')
% plot(dz1, 'b--')


% figure
% hold on
% plot(mdx2, 'r')
% plot(mdy2, 'g')
% plot(mdz2, 'b')
% plot(dx2, 'r--')
% plot(dy2, 'g--')
% plot(dz2, 'b--')


