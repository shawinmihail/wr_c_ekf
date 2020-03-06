% clc
% clear
close all

%% est
path = 'est_state_log.csv';
data = csvread(path);

c_x_est = data(:,1);
c_y_est = data(:,2);
c_z_est = data(:,3);

c_vx_est = data(:,4);
c_vy_est = data(:,5);
c_vz_est = data(:,6);

c_ax_est = data(:,7);
c_ay_est = data(:,8);
c_az_est = data(:,9);

c_qw_est = data(:,10);
c_qx_est = data(:,11);
c_qy_est = data(:,12);
c_qz_est = data(:,13);

c_wx_est = data(:,14);
c_wy_est = data(:,15);
c_wz_est = data(:,16);

% figure
% hold on
% plot(c_x_est, 'r')
% plot(c_y_est, 'g')
% plot(c_z_est, 'b')

run('parse_sim_data');


figure
hold on
plot(x_est, 'r--')
plot(y_est, 'g--')
plot(z_est, 'b--')

plot(c_x_est, 'r')
plot(c_y_est, 'g')
plot(c_z_est, 'b')

figure
hold on
plot(vx_est, 'r--')
plot(vy_est, 'g--')
plot(vz_est, 'b--')

plot(c_vx_est, 'r')
plot(c_vy_est, 'g')
plot(c_vz_est, 'b')

figure
hold on
plot(qx_est, 'r--')
plot(qy_est, 'g--')
plot(qz_est, 'b--')

plot(c_qx_est, 'r')
plot(c_qy_est, 'g')
plot(c_qz_est, 'b')



