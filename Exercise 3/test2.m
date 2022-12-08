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
        
        h_bar_2 = zeros(N,1);
        for n = 1 : N
            for l = 0 : L-1
                h_bar_2(n) = h_bar_2(n) + h(l+1)*exp( -(1i*2*pi*l*n)/N );
            end
        end
        h_bar_2 = (1/sqrt(N))*h_bar_2;
        
        figure;
        scatter(real(h_bar),imag(h_bar));
        hold on 
        scatter(real(h_bar_2),imag(h_bar_2));
        
        % 4-QAM data block
        d = sign(-1+2*rand(N,1)) + 1i*sign(-1+2*rand(N,1));
        d_bar = (1/sqrt(N))*fft(d,N);
        
        d_bar_2 = zeros(N,1);
        for n = 1:N
            for m = 0:N-1
                d_bar_2(n) = d_bar_2(n) + d(m+1)*exp( -(1i*2*pi*m*n)/N );
            end
        end
        d_bar_2 = (1/sqrt(N))*d_bar_2;
        
        figure;
        scatter(real(d_bar),imag(d_bar));
        hold on 
        scatter(real(d_bar_2),imag(d_bar_2));
        
        % Channel input
        x = [d(N-L+2 : N) ; d];
        
        % Channel noise
        w = ( randn(N,1) + 1i*randn(N,1) )*sqrt(N_0/2);
        w_bar = (1/sqrt(N))*fft(w,N);
        
        % Received Signal
        y_1 = zeros(N,1);
        for m = 1 : N
            for l = 1 : L
                y_1(m) = y_1(m) + h(l)*x( m+L-1-(l-1) );
            end
        end
        
        y_2 = zeros(N+L-1,1);
        for m = L : N+L-1
            for l = 0 : L-1
                y_2(m) = y_2(m) + h(l+1)*x(m-l);
            end
        end
        y_2 = y_2(L : end);
        
%         y_3 = zeros(N,1);
%         for m = 1 : N
%             y_3(m) = h(1)*d(m);
%             for l = 1 : L-1
%                 y_3(m) = y_3(m) + h(l+1)*d(N-l);
%             end
%         end
        
        figure;
        scatter(real(y_1),imag(y_1));
        hold on
        scatter(real(y_2),imag(y_2));
%         scatter(real(y_3)+0.2,imag(y_3)+0.2);
        
        y_bar = sqrt(N)*h_bar.*d_bar;
        y_bar_2 = sqrt(N)*h_bar_2.*d_bar_2;
        
        figure;
        scatter(real(y_1),imag(y_1));
        hold on 
        scatter(real(y_2),imag(y_2));
        
        figure;
        scatter(real(y_bar),imag(y_bar));
        hold on 
        scatter(real(y_bar_2),imag(y_bar_2));

%         y_2 = y_2 + w;
%         y_bar = y_bar + w_bar;
        
%         figure;
%         scatter(real(y_ifft_0),imag(y_ifft_0));
%         hold on 
%         scatter(real(y_bar),imag(y_bar));
%         scatter(real(y_bar_2),imag(y_bar_2));
        
        y_ifft = sqrt(N)*ifft(y_bar,N);
        y_ifft_2 = sqrt(N)*ifft(y_bar_2,N);
        figure;
        scatter(real(y_2),imag(y_2));
        hold on 
        scatter(real(y_ifft),imag(y_ifft));
        scatter(real(y_ifft_2),imag(y_ifft_2));
        
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
