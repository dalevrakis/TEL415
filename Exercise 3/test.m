clc; close all; clear;

K = 2; % User Number
M = 100; % 
N = 64; % Sequence lengh
L = 4; % Responce Length
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
        h_bar = (1/sqrt(N))*fft(h,N);
        
        % 4-QAM data block
        d = sign(-1+2*rand(N,1)) + 1i*sign(-1+2*rand(N,1));
        d_bar = (1/sqrt(N))*fft(d,N);
        
        % Channel input
        x = [d(N-L+2 : N) ; d];
        
        % Channel noise
        w = ( randn(N,1) + 1i*randn(N,1) )*sqrt(N_0/2);
        w_bar = (1/sqrt(N))*fft(w,N);
        
        % Received Signal
%         y = zeros(N,1);
%         for m = 1 : N
%             for l = 1 : L
%                 y(m) = y(m) + h(l)*x( m+L-1-(l-1) );
%             end
%         end
        
        y = zeros(N+L-1,1);
        for m = L : N+L-1
            for l = 0 : L-1
                y(m) = y(m) + h(l+1)*x(m-l);
            end
        end
        y = y(L : end);
        y_bar = sqrt(N)*h_bar.*d_bar;
        
        figure;
        scatter(real(y),imag(y));
        hold on 
        scatter(real(y_bar),imag(y_bar));
        
        y = y + w;
        y_bar = y_bar + w_bar;
        
        figure;
        scatter(real(y),imag(y));
        hold on 
        scatter(real(y_bar),imag(y_bar));
        
        y_ifft = (1/sqrt(N))*ifft(y_bar,N);
%         y_dist = zeros(4,N);
%         y_dist(1,:) = abs(y-X(1)).^2;
%         y_dist(2,:) = abs(y-X(2)).^2;
%         y_dist(3,:) = abs(y-X(3)).^2;
%         y_dist(4,:) = abs(y-X(4)).^2;
% 
%         decision_matrix = zeros(N,1);
%         for k = 1:N
%             min_dist = find( y_dist(:,k) == min(y_dist(:,k)) );
%             decision_matrix(k) = X(min_dist(1));
%         end
% 
%         error_matrix = d-decision_matrix;
%         error_indices = find(error_matrix);
%         
%         err_bits = err_bits + sum( abs( (error_matrix).^2 )/4 );
    end
    
    BER_snr(snr) = err_bits/(N*Packets);
end
toc;
