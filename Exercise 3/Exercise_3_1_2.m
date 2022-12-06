clc; close all; clear;

K = 1; % User Number
M = 300; % Sequence lengh
N = 64; % Code Length
L = 3; % Responce Length
Packets = 1000;
X = [1, -1];

snr_db = 0: 2 : 20;
BER_snr = zeros(size(snr_db,2),1);
for snr = 1 : size(snr_db,2)
    %AWGN regarding SNR
    N_0 = 1/( 10^(snr_db(snr)/10) );
    
    err_bits = 0;
    for pac = 1:Packets
        %User code
        c = (1/sqrt(N))*sign(randn(N,K));

        %User Symbol Sequence
        s = sign(-1+2*rand(M,K));

        c_l = zeros(N+L-1,L);
        for l = 1:L
            c_l(l:N+l-1,l) = c; 
        end

        %Channel response
        h = (randn(1,L) + 1i*randn(1,L))*sqrt(1/(2*L));
        [~,l_max] = max(h);
        norm_h = norm(h,2);
        
        r = zeros(M,1);
        decision_matrix = zeros(M,1);
        error_matrix = zeros(M,1);
        for m = 1 : M
            x_m = s(m)*c_l;
            y = zeros(N+L-1,1);
            for l = 1 : L
                y = y + h(l)*x_m(:,l);
            end
            y = y + ( randn(N+L-1,1) + 1i*randn(N+L-1,1) )*sqrt(N_0/2);

            % Rake Receiver
%             for l = 1 : L
%                 r(m) = r(m) + sum( (conj(h(l))/norm_h)*c_l(:,l).*y );
%             end
            r(m) = sum( (conj(h(l_max))/norm_h)*c_l(:,l_max).*y );
            decision_matrix(m) = sign(real(r(m)));

            error_matrix(m) = s(m)-decision_matrix(m);
        end
        err_bits = err_bits+sum( abs( (error_matrix).^2 )/4 );
        
    end
    BER_snr(snr) = err_bits/(M*Packets);
%     figure;
%     scatter(real(r),imag(r));
end
figure;
semilogy(snr_db,BER_snr);
hold on
fi = 0:0.1:20;
semilogy(fi,1./(fi));
semilogy(fi,1./(fi.^2));
semilogy(fi,1./(fi.^3));
