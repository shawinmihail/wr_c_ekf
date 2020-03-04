clc
clear
close all

%% act
path = 'act_state_log.csv';
data = csvread(path);

x_act = data(:,1);
y_act = data(:,2);
z_act = data(:,3);

vx_act = data(:,4);
vy_act = data(:,5);
vz_act = data(:,6);

ax_act = data(:,7);
ay_act = data(:,8);
az_act = data(:,9);

qx_act = data(:,10);
qy_act = data(:,11);
qz_act = data(:,12);

wx_act = data(:,13);
wy_act = data(:,14);
wz_act = data(:,15);

%% mes
path = 'mes_state_log.csv';
data = csvread(path);

x_mes = data(:,1);
y_mes = data(:,2);
z_mes = data(:,3);

vx_mes = data(:,4);
vy_mes = data(:,5);
vz_mes = data(:,6);

ax_mes = data(:,7);
ay_mes = data(:,8);
az_mes = data(:,9);

qx_mes = data(:,10);
qy_mes = data(:,11);
qz_mes = data(:,12);

wx_mes = data(:,13);
wy_mes = data(:,14);
wz_mes = data(:,15);

%% est
path = 'est_state_log.csv';
data = csvread(path);

x_est = data(:,1);
y_est = data(:,2);
z_est = data(:,3);

vx_est = data(:,4);
vy_est = data(:,5);
vz_est = data(:,6);

ax_est = 0;
ay_est = 0;
az_est = 0;

qx_est = data(:,7);
qy_est = data(:,8);
qz_est = data(:,9);

wx_est = 0;
wy_est = 0;
wz_est = 0;


figure
hold on
plot(qx_est, 'r')
plot(qy_est, 'g')

figure
hold on
plot(x_mes, y_mes,  'y')
plot(x_est, y_est, 'b')
% plot(x_act, y_act, 'r')


% plot(vx_act, 'k')

