y1 = readtable('y1.txt', 'ReadVariableNames', false);
y2 = readtable('y2.txt', 'ReadVariableNames', false);
y1_fin = table2array(y1);
y2_fin = table2array(y2);
t = 0:0.001:1-0.001;
ni = 2.0;
y1_tacno = sin(t) + ni .* cos(t) .* t;
y2_tacno = -cos(t);
plot(t,y1_fin);
hold on;
plot(t,y1_tacno);
grid on;
legend('Postinuta vrijednost','Egzaktna vrijednost');
title('Vrijednost $y_1$');
xlabel('$t(s)$');
ylabel('$y_1$');
figure(2)
plot(t,y2_fin);
hold on;
plot(t,y2_tacno);
grid on;
legend('Postinuta vrijednost','Egzaktna vrijednost');
title('Vrijednost $y_2$');
xlabel('$t(s)$');
ylabel('$y_2$');

figure(3)
plot(t,y1_fin' - y1_tacno);
grid on;
title('Vrijednost greske po $y_1$');
xlabel('$t(s)$');
ylabel('$e_{y_1}$');

figure(4)
plot(t,y2_fin' - y2_tacno);
grid on;
title('Vrijednost greske po $y_2$');
xlabel('$t(s)$');
ylabel('$e_{y_2}$');