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
        h = zeros(M,N);
        for div = 1 : M
            h(div,:) = (randn(1,N) + 1i*randn(1,N))*sqrt(1/2);
        end
        
        %% 2
        s = 1 - (randi(2,[1,N])-1)*2 + 1i*(  1 - (randi(2,[1,N])-1)*2 );

        n = zeros(M,N);
        for div = 1:M
            n(div,:) = ( randn(1,N) + 1i*randn(1,N) )*sqrt(N_0/2);
        end

        r = zeros(M,N);
        for div = 1:M
            r(div,:) = h(div,:).*s + n(div,:);
        end

        %% 4

        R = zeros(1,N);
        for i = 1:N
            R(i) = (h(:,i)'/norm(h(:,i),2))*r(:,i);
        end

        y_dist = zeros(4,N);
        y_dist(1,:) = abs(R-X(1)).^2;
        y_dist(2,:) = abs(R-X(2)).^2;
        y_dist(3,:) = abs(R-X(3)).^2;
        y_dist(4,:) = abs(R-X(4)).^2;

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
% fig1=figure;
% semilogy(snr_db,snr_BER)
% xlabel('$SNR_{db}$','Interpreter','latex');
% ylabel('BER','Interpreter','latex');
% % title('Real part visualization of h');
% saveas(fig1,'fig1.png')

fig4=figure;
semilogy(snr_db,snr_BER)
hold on
semilogy(snr_db,snr_Pe)

for snr = 1 : snr_size
    err_bits = 0;
    N_0 = 2/10^(snr_db(snr)/10);
    for packet = 1 : K
        %% 1
        h = zeros(M,N);
        for div = 1 : M
            h(div,:) = (randn(1,N) + 1i*randn(1,N))*sqrt(1/2);
        end
        
        %% 2
        s = 1 - (randi(2,[1,N])-1)*2 + 1i*(  1 - (randi(2,[1,N])-1)*2 );
        
%         n = zeros(M,N);
%         for div = 1:M
%             n(div,:) = ( randn(1,N) + 1i*randn(1,N) )*sqrt(N_0/2);
%         end
        n = ( randn(1,N) + 1i*randn(1,N) )*sqrt(N_0/2);
        
%         r = zeros(M,N);
%         for div = 1:M
%             r(div,:) = h(div,:).*s + n(div,:);
%         end
        
        r = zeros(1,N);
        for i = 1:N
            r(i) = norm(h(:,i),2)*s(i)+n(i);
        end
        %% 4

%         R = zeros(1,N);
%         for i = 1:N
%             R(i) = (h(:,i)'/norm(h(:,i),2))*r(:,i);
%         end
% 
%         y_dist = zeros(4,N);
%         y_dist(1,:) = abs(R-X(1)).^2;
%         y_dist(2,:) = abs(R-X(2)).^2;
%         y_dist(3,:) = abs(R-X(3)).^2;
%         y_dist(4,:) = abs(R-X(4)).^2;
% 
%         decision_matrix = zeros(1,N);
%         for k = 1:N
%             min_dist = find( y_dist(:,k) == min(y_dist(:,k)) );
%             decision_matrix(k) = X(min_dist(1));
%         end

        y_dist = zeros(4,N);
        y_dist(1,:) = abs(r-X(1)).^2;
        y_dist(2,:) = abs(r-X(2)).^2;
        y_dist(3,:) = abs(r-X(3)).^2;
        y_dist(4,:) = abs(r-X(4)).^2;

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

semilogy(snr_db,snr_BER)

xlabel('$SNR_{db}$','Interpreter','latex');
ylabel('BER','Interpreter','latex');
legend({'MRC BER','Theoretical BER','TB BER'},'Interpreter','latex');
saveas(fig4,'fig4.png')
legend show