clc; close all; clear;

desired_snr_db = 3;
Rh0 = 1*10^(desired_snr_db/10)/2;
b = 0.999;
N = max(1000,ceil(1/(1-b))*100);
%% 1
h = zeros(N,1);
h(1) = 1 + 1i*0;

sigma_e2 = Rh0*(1-b^2);

e = (randn(N,1) + 1i*randn(N,1))*sqrt(sigma_e2/2);


for k = 2:N
    h(k) = b*real(h(k-1)) + real(e(k)) + 1i*( b*imag(h(k-1)) + imag(e(k)) );
end

figure;
plot(1:N,real(h));
xlabel('\bf{k}','Interpreter','latex');
ylabel('\bf{R[h]}','Interpreter','latex');
title('Real part visualization of h');

figure;
plot(1:N,imag(h));
xlabel('\bf{k}','Interpreter','latex');
ylabel('\bf{I[h]}','Interpreter','latex');
title('Imaginary part visualization of h');
%% 2
s = 1 - (randi(2,[N,1])-1)*2 + 1i*(  1 - (randi(2,[N,1])-1)*2 );

%% 3
n = ( randn(N,1) + 1i*randn(N,1) )*sqrt(1/2);

r = h.*s + n;

snr = mean(abs(h.*s).^2)/mean(abs(n).^2);
snr_db  = 10*log10( snr );

%% 4
y = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i ];

y_dist = zeros(N,4);
y_dist(:,1) = abs(r-h.*y(1)).^2;
y_dist(:,2) = abs(r-h.*y(2)).^2;
y_dist(:,3) = abs(r-h.*y(3)).^2;
y_dist(:,4) = abs(r-h.*y(4)).^2;

decision_matrix = zeros(N,1);
for k = 1:N
    min_dist = find( y_dist(k,:) == min(y_dist(k,:)) );
    decision_matrix(k) = y(min_dist(1));
end

error_matrix = s-decision_matrix;
error_indices = find(error_matrix);
exper_P_err = size( error_indices,1 )/N;
theor_P_err = 1/(2*snr);

%% 5

figure;
subplot(211);
semilogy(1:N, abs(h));
ylabel('$\bf{|h|}$','Interpreter','latex');
title('Comparison of magnitude of h and errors');

subplot(212);
stem(1:N,abs(error_matrix).^2,"|");
xlabel('\bf{k}','Interpreter','latex');
ylabel('$|error|^2$','Interpreter','latex');