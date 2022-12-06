clc; close all; clear;

N = 100;
M = 2;
X = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i ];
K = 1000;

%% 3
snr_db = 0 : 2 : 20;
snr_size = size(snr_db,2);
snr_BER = zeros(snr_size,1);
snr_Pe = zeros(snr_size,1);

for snr = 1 : snr_size
    err_bits = 0;
    N_0 = 2/10^(snr_db(snr)/10);
    for packet = 1 : K
        %% 1
        h_1 = (randn(1,N/2) + 1i*randn(1,N/2))*sqrt(1/2);
        h_2 = (randn(1,N/2) + 1i*randn(1,N/2))*sqrt(1/2);
        
        %% 2
        s = 1 - (randi(2,[1,N])-1)*2 + 1i*(  1 - (randi(2,[1,N])-1)*2 );
        

        n_1 = ( randn(1,N/2) + 1i*randn(1,N/2) )*sqrt(N_0/2);
        n_2 = ( randn(1,N/2) + 1i*randn(1,N/2) )*sqrt(N_0/2);
        
        r_1 = zeros(1,N/2);
        r_2 = zeros(1,N/2);
        for i = 1:N/2
            r_1(i) = h_1(i)*s(2*i-1) + h_2(i)*s(2*i) + n_1(i);
            r_2(i) = h_2(i)*conj(s(2*i-1)) - h_1(i)*conj(s(2*i)) + n_2(i);
        end
        
        %% 4
          R_1 = zeros(1,N/2);
          R_2 = zeros(1,N/2);
          for i = 1:N/2
            h_b = [h_1(i) h_2(i)]';
            r = [r_1(i) ; conj(r_2(i))];
            
            H = (1/norm(h_b,2))*[h_1(i) h_2(i) ; conj(h_2(i)) -conj(h_1(i)) ];
            H_H = H';
            R_1(i) = H_H(1,:)*r;
            R_2(i) = H_H(2,:)*r;
          end
          
          y_dist = zeros(4,N);
          for i = 1:N/2
            h_b = [h_1(i) h_2(i)]';
            
            y_dist(1,2*i-1) = abs(R_1(i)-norm(h_b,2)*X(1)).^2;
            y_dist(2,2*i-1) = abs(R_1(i)-norm(h_b,2)*X(2)).^2;
            y_dist(3,2*i-1) = abs(R_1(i)-norm(h_b,2)*X(3)).^2;
            y_dist(4,2*i-1) = abs(R_1(i)-norm(h_b,2)*X(4)).^2;
            
            y_dist(1,2*i) = abs(R_2(i)-norm(h_b,2)*X(1)).^2;
            y_dist(2,2*i) = abs(R_2(i)-norm(h_b,2)*X(2)).^2;
            y_dist(3,2*i) = abs(R_2(i)-norm(h_b,2)*X(3)).^2;
            y_dist(4,2*i) = abs(R_2(i)-norm(h_b,2)*X(4)).^2;
          end

        decision_matrix = zeros(1,N);
        for k = 1:N
            min_dist = find( y_dist(:,k) == min(y_dist(:,k)) );
            decision_matrix(k) = X(min_dist(1));
        end


        error_matrix = s-decision_matrix;
        error_indices = find(error_matrix);
        
        err_bits = err_bits + sum( abs( (error_matrix).^2 )/4 );
    end
    snr_BER(snr) = err_bits/(N*K);
    
    snr_Pe(snr) = nchoosek(2*M-1,M)*( 1/( (4^M) * snr_db(snr)^M ) );
end
fig5=figure;
semilogy(snr_db,snr_BER)

xlabel('$SNR_{db}$','Interpreter','latex');
ylabel('BER','Interpreter','latex');
legend({'Alamouti BER'},'Interpreter','latex');
saveas(fig5,'fig5.png')
legend show
% hold on
% semilogy(snr_db,snr_Pe)