
korak = 0.05;
vrijeme=10;
y1 = readtable('y1.txt', 'ReadVariableNames', false);
y2 = readtable('y2.txt', 'ReadVariableNames', false);
w1 = readtable('w1.txt', 'ReadVariableNames', false);
w2 = readtable('w2.txt', 'ReadVariableNames', false);
w3 = readtable('w3.txt', 'ReadVariableNames', false);
delta1 = readtable('delta1.txt', 'ReadVariableNames', false);
Vr1 = readtable('Vr1.txt', 'ReadVariableNames', false);
psv1 = readtable('psv1.txt', 'ReadVariableNames', false);
edp1 = readtable('Edp1.txt', 'ReadVariableNames', false);
eqp1= readtable('Eqp1.txt', 'ReadVariableNames', false)


y1_fin = table2array(y1);
y2_fin = table2array(y2);
w1_fin = table2array(w1);
w2_fin = table2array(w2);
w3_fin = table2array(w3);
delta1_fin = table2array(delta1);
vr1_fin = table2array(Vr1);
psv1_fin = table2array(psv1);
edp1_fin = table2array(edp1);
eqp1_fin = table2array(eqp1);


w1_fin = w1_fin / (2 * pi);
w2_fin = w2_fin / (2 * pi);
w3_fin = w3_fin / (2 * pi);


t = 0:korak:vrijeme-korak;
figure(1)
plot(t,w1_fin);
grid on;
xlabel('$t(s)$');
ylabel('$\omega_1$');

figure(2)
plot(t,w2_fin);
grid on;
xlabel('$t(s)$');
ylabel('$\omega_2$');


figure(3)
plot(t,w3_fin);
grid on;
xlabel('$t(s)$');
ylabel('$\omega_3$');


figure(4)

plot(t,delta1_fin);
xlim([0 10]);
ylim([-0.4 0.8]);
grid on;
xlabel('$t(s)$');
ylabel('$\delta_1$');


figure(5)
plot(t,vr1_fin);
grid on;
xlabel('$t(s)$');
ylabel('$Vr_1$');

figure(6)
plot(t,psv1_fin);
grid on;
xlabel('$t(s)$');
ylabel('$PSV_1$');

figure(7)
plot(t,edp1_fin);
grid on;
xlabel('$t(s)$');
ylabel('$Edp_1$');



figure(8)
plot(t,eqp1_fin);
grid on;
xlabel('$t(s)$');
ylabel('$Eqp_1$');
