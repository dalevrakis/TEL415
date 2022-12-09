clc; close all; clear;

K = 5; % User Number
M = 100; % Sequence lengh
N = 128; % Code Length
L = 3; % Responce Length
Packets = 1000;
X = [1, -1];

snr_db = 0: 2 : 16;
BER_snr = zeros(size(snr_db,2),1);
tic;
for snr = 1 : size(snr_db,2)
    %AWGN regarding SNR
    N_0 = 1/( 10^(snr_db(snr)/10) );
    
    err_bits = 0;
    for pac = 1:Packets
        %User code
        c = (1/sqrt(N))*sign(randn(N,K));

        %User Symbol Sequence
        s = sign(-1+2*rand(M,K));
        
        c_l = zeros(K,N+L-1,L);
        for k = 1:K
            for l = 1:L
                c_l(k,l:N+l-1,l) = c(:,k); 
            end
        end

        %Channel response
        h = (randn(K,L) + 1i*randn(K,L))*sqrt(1/(2*L));
        norm_h1 = norm(h(1,:),2);
        
        r = zeros(M,1);
        decision_matrix = zeros(M,1);
        error_matrix = zeros(M,1);
        for m = 1 : M
            x_m = zeros(K,N+L-1,L);
            for k = 1:K
                x_m(k,:,:) = s(m,k)*c_l(k,:,:);
            end
          
            y = zeros(N+L-1,1);
            for k = 1:K
                for l = 1 : L
                    y = y + h(k,l)*x_m(k,:,l)';
                end
            end
            y = y + ( randn(N+L-1,1) + 1i*randn(N+L-1,1) )*sqrt(N_0/2);

            % Rake Receiver
            for l = 1 : L
                r(m) = r(m) + sum( (conj(h(1,l))/norm_h1)*c_l(1,:,l)'.*y );
            end

            decision_matrix(m) = sign(real(r(m)));

            error_matrix(m) = s(m,1)-decision_matrix(m);
        end
        err_bits = err_bits+sum( abs( (error_matrix).^2 )/4 );
        
    end
    BER_snr(snr) = err_bits/(M*Packets);
%     figure;
%     scatter(real(r),imag(r));
end
toc;
fig10=figure;
semilogy(snr_db,BER_snr);
hold on
fi = 0:0.1:16;
semilogy(fi,1./(fi));
semilogy(fi,1./(fi.^2));
semilogy(fi,1./(fi.^3));
xlabel('$SNR_{db}$','Interpreter','latex');
ylabel('BER','Interpreter','latex');
legend({'CDMA','$\frac{1}{SNR}$','$\frac{1}{SNR^2}$','$\frac{1}{SNR^3}$'},'Interpreter','latex');
saveas(fig10,'fig10.png')
legend show