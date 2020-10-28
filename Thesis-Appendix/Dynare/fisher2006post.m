dynare fisher2006.mod;

a_c_n = [0; cumsum(a_ea(1:end-1))];
v_c_n = [0; cumsum(v_ea(1:end-1))];
a_c_i = [0; cumsum(a_ev(1:end-1))];
v_c_i = [0; cumsum(v_ev(1:end-1))];

%a_c_n = cumsum(a_ea(1:end));
%v_c_n = cumsum(v_ea(1:end));
%a_c_i = cumsum(a_ev(1:end));
%v_c_i = cumsum(v_ev(1:end));


y_c_n = y_ea + a_c_n/(1-alpha) + v_c_n*alpha/(1-alpha);
y_c_i = y_ev + a_c_i/(1-alpha) + v_c_i*alpha/(1-alpha);

h_c_n = h_ea;
h_c_i = h_ev;

c_c_n = c_ea + a_c_n/(1-alpha) + v_c_n*alpha/(1-alpha);
c_c_i = c_ev + a_c_i/(1-alpha) + v_c_i*alpha/(1-alpha);

x_c_n = x_ea + a_c_n/(1-alpha) + v_c_n*alpha/(1-alpha);
x_c_i = x_ev + a_c_i/(1-alpha) + v_c_i*alpha/(1-alpha);

k_c_n = k_ea + a_c_n/(1-alpha) + v_c_n/(1-alpha);
k_c_i = k_ev + a_c_i/(1-alpha) + v_c_i/(1-alpha);


% 
% % N-shock
% figure(11);
% suptitle('N-shock');
% 
% subplot(321)
% sr = -v_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Pt'); hold off;
% subplot(323)
% sr = h_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Ht'); hold off;
% subplot(325)
% sr = c_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Ct'); hold off;
% subplot(322)
% sr = y_c_n - h_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Yt/Ht'); hold off;
% subplot(324)
% sr = y_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:');  title('Yt'); hold off;
% subplot(326)
% sr = v_c_n + x_c_n;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:');  title('VtXt'); hold off;
% 
% figure(12);
% suptitle('I-shock');
% 
% subplot(321)
% sr = -v_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Pt'); hold off;
% subplot(323)
% sr = h_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Ht'); hold off;
% subplot(325)
% sr = c_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Ct'); hold off;
% subplot(322)
% sr = y_c_i - h_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:'); title('Yt/Ht'); hold off;
% subplot(324)
% sr = y_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:');  title('Yt'); hold off;
% subplot(326)
% sr = v_c_i + x_c_i;
% hold on; plot(0:H, [0; sr]); plot([0 H], [0 0], 'k:');  title('VtXt'); hold off;

% TODO plot both
figure(13);

subplot(321)
sr1 = -v_c_i; sr2 = -v_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:'); title('Pt'); hold off;
subplot(323)
sr1 = h_c_i; sr2 = h_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:'); title('Ht'); hold off;
subplot(325)
sr1 = c_c_i; sr2 = c_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:');  title('Ct'); hold off;
subplot(322)
sr1 = y_c_i - h_c_i; sr2 = y_c_n - h_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:'); title('Yt/Ht'); hold off;
subplot(324)
sr1 = y_c_i; sr2 = y_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:');  title('Yt'); hold off;
subplot(326)
sr1 = v_c_i + x_c_i; sr2 = v_c_n + x_c_n;
hold on; plot(0:H, [0; sr1], 'b--'); plot(0:H, [0; sr2], 'r-'); plot([0 H], [0 0], 'k:'); title('VtXt'); hold off;
legend('I','N')