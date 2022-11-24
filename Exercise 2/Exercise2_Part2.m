clc; close all; clear;

N = 100;
i = 2;
j = 2;
X_s = [1 + 1i, -1 + 1i, 1 - 1i, -1 - 1i ];
K = 1000;

%% 3
snr_db = 0 : 2 : 20;
snr_size = size(snr_db,2);
snr_BER_ml = zeros(snr_size,1);
snr_BER_decol = zeros(snr_size,1);
snr_Pe = zeros(snr_size,1);

for snr = 1 : snr_size
    err_bits_ml = 0;
    err_bits_decol = 0;
    N_0 = 2/10^(snr_db(snr)/10);
    for packet = 1 : K
        %% 1
        h = zeros(2,2,N);
        for div_i = 1 : i
            for div_j = 1 : j
                h(div_i,div_j,:) = (randn(1,N) + 1i*randn(1,N))*sqrt(1/2);
            end 
        end
        
        %% 2
        X = zeros(i,N);
        for div_i = 1:i
            X(div_i,:) = 1 - (randi(2,[1,N])-1)*2 + 1i*(  1 - (randi(2,[1,N])-1)*2 );
        end
        
        %% 3
        W = zeros(i,N);
        for div_i = 1:i
            W(div_i,:) = ( randn(1,N) + 1i*randn(1,N) )*sqrt(N_0/2);
        end
        
        Y = zeros(div_i,N);
        
%         for div_i = 1:i
%             Y(div_i,:) = squeeze(h(div_i,1,:))'.*X(1,:) + squeeze(h(div_i,2,:))'.*X(2,:) + W(div_i,:);
%         end
        
        for t = 1:N
            Y(:,t) = h(:,:,t)*X(:,t) + W(:,t);
        end
        
        %% 4 a
        combinations_ml = size(X_s,2)*size(X_s,2);
        X_comb = zeros(2,combinations_ml);
        for comb_1 = 1 : size(X_s,2)
            for comb_2 = 1 : size(X_s,2)
                X_comb(1,(comb_1-1)*size(X_s,2) + comb_2) = X_s(comb_1);
                X_comb(2,(comb_1-1)*size(X_s,2) + comb_2) = X_s(comb_2);
            end
        end
        
        y_dist = zeros(combinations_ml,N);
        for t = 1 : N
            for comb = 1 : combinations_ml
                y_dist(comb,t) = norm(Y(:,t)-squeeze(h(:,:,t))*X_comb(:,comb),2).^2;
            end
        end
       

        decision_matrix_ml = zeros(2,N);
        for k = 1:N
            min_dist = find( y_dist(:,k) == min(y_dist(:,k)) );
            decision_matrix_ml(:,k) = X_comb(:,min_dist(1));
        end


        error_matrix_ml = X-decision_matrix_ml;
        error_indices_ml = find(error_matrix_ml);
        
        err_bits_ml = err_bits_ml + sum( abs( (error_matrix_ml(1,:)).^2 )/4 );
        err_bits_ml = err_bits_ml + sum( abs( (error_matrix_ml(2,:)).^2 )/4 );
        
        
        %% 4 b
        X_tilde = zeros(2,N);
        for t = 1 : N
            X_tilde(:,t) = h(:,:,t)\Y(:,t);
        end
        
        combinations_decol = size(X_s,2);
        y_dist_1 = zeros(combinations_decol,N);
        for t = 1 : N
            for comb = 1 : combinations_decol
                y_dist_1(comb,t) = norm(X_tilde(1,t)-X_s(comb),2).^2;
            end
        end
       
        y_dist_2 = zeros(combinations_decol,N);
        for t = 1 : N
            for comb = 1 : combinations_decol
                y_dist_2(comb,t) = norm(X_tilde(2,t)-X_s(comb),2).^2;
            end
        end
        
        decision_matrix_decol = zeros(2,N);
        for k = 1:N
            min_dist = find( y_dist_1(:,k) == min(y_dist_1(:,k)) );
            decision_matrix_decol(1,k) = X_s(min_dist(1));
            
            min_dist = find( y_dist_2(:,k) == min(y_dist_2(:,k)) );
            decision_matrix_decol(2,k) = X_s(min_dist(1));
        end


        error_matrix_decol = X-decision_matrix_decol;
        error_indices_decol = find(error_matrix_decol);
        
        err_bits_decol = err_bits_decol + sum( abs( (error_matrix_decol(1,:)).^2 )/4 );
        err_bits_decol = err_bits_decol + sum( abs( (error_matrix_decol(2,:)).^2 )/4 );
    end
    snr_BER_ml(snr) = err_bits_ml/(N*K);
    
    snr_BER_decol(snr) = err_bits_decol/(N*K);
end
figure;
semilogy(snr_db,snr_BER_ml)
hold on
semilogy(snr_db,snr_BER_decol)