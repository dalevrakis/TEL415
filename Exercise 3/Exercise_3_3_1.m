clc; close all; clear;

N = 128; % Sequence lengh
L = 8; % Responce Length
Packets = 1;
X = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i ];


snr_db = 15;
BER_snr = zeros(size(snr_db,2),1);
tic;
for snr = 1 : size(snr_db,2)
    %AWGN regarding SNR
    N_0 = 1/( 10^(snr_db(snr)/10) );
    
    err_bits = 0;
    for pac = 1:Packets
        % Channel response
        h = (randn(L,1) + 1i*randn(L,1))*sqrt(1/(2*L));
        h_tilde = (1/sqrt(N))*fft(h,N);
        
        % 4-QAM data block
        d = sign(-1+2*rand(N,1)) + 1i*sign(-1+2*rand(N,1));
        d_tilde = (1/sqrt(N))*fft(d,N);

        % Channel input
        x = [d(N-L+2 : N) ; d];
        
        % Channel noise
        w = ( randn(N,1) + 1i*randn(N,1) )*sqrt(N_0/2);
        w_tilde = (1/sqrt(N))*fft(w,N);
        
        % Received Signal
        y = zeros(N+L-1,1);
        for m = L : N+L-1
            for l = 0 : L-1
                y(m) = y(m) + h(l+1)*x(m-l);
            end
        end
        y = y(L : end);
        
        y_tilde = sqrt(N)*h_tilde.*d_tilde;
        
        %N-Parallel - Frequency Selective Comparison
        y_freq = sqrt(N)*ifft(y_tilde,N);
        
        figure;
        scatter(real(y),imag(y),'o');
        hold on
        scatter(real(y_freq),imag(y_freq),'x');
        
        y = y + w;
        y_tilde = y_tilde + w_tilde;
        y_freq_noise = sqrt(N)*ifft(y_tilde);
        
        figure;
        scatter(real(y),imag(y),'o');
        hold on
        scatter(real(y_freq_noise),imag(y_freq_noise),'x');
        
        
        y_dec = y_tilde./h_tilde;
        d_rec = ifft(y_dec,N);
        figure;
        scatter(real(d_rec),imag(d_rec));
        hold on
        scatter(real(d),imag(d));
        
        %Decision
        d_dist = zeros(4,N);
        d_dist(1,:) = abs(d_rec-X(1)).^2;
        d_dist(2,:) = abs(d_rec-X(2)).^2;
        d_dist(3,:) = abs(d_rec-X(3)).^2;
        d_dist(4,:) = abs(d_rec-X(4)).^2;

        decision_matrix = zeros(N,1);
        for k = 1:N
            min_dist = find( d_dist(:,k) == min(d_dist(:,k)) );
            decision_matrix(k) = X(min_dist(1));
        end

        error_matrix = d-decision_matrix;
        error_indices = find(error_matrix);
        
        err_bits = err_bits + sum( abs( (error_matrix).^2 )/4 );
        
        figure;
        plot(1:N,abs(h_tilde))
        hold on
        plot(1:N,abs(w_tilde))
        stem(error_indices,abs(h_tilde(error_indices)))
        
        figure;
        plot(1:N,abs(h_tilde))
        hold on
        plot(1:N,abs(w_tilde))
        stem(error_indices,abs(w_tilde(error_indices)))
    end
    
    BER_snr(snr) = err_bits/(N*Packets);
end
toc;
