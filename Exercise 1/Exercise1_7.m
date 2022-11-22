clc; close all; clear;

desired_snr_db = 0 : 2 : 30;
snr_err_bits = zeros(size(desired_snr_db,1),1);
snr_BER = zeros(size(desired_snr_db,1),1);
packet_num = 5000;
b = 0.4;
N = max(200,ceil(1/(1-b))*10);
    
for snr_i = 1 : size(desired_snr_db,2)
    Rh0 = 1*10^(desired_snr_db(snr_i)/10)/2;
    sigma_e2 = Rh0*(1-b^2);
    err_bits = 0;
    
%     h = zeros(N,1);
%     h(1) = 1 + 1i*0;
% 
%     e = (randn(N,1) + 1i*randn(N,1))*sqrt(sigma_e2/2);
% 
%     for k = 2:N
%         h(k) = b*real(h(k-1)) + real(e(k)) + 1i*( b*imag(h(k-1)) + imag(e(k)) );
%     end
    for P = 1 : packet_num
        %% 1
        h = zeros(N,1);
        h(1) = 1 + 1i*0;

        e = (randn(N,1) + 1i*randn(N,1))*sqrt(sigma_e2/2);

        for k = 2:N
            h(k) = b*real(h(k-1)) + real(e(k)) + 1i*( b*imag(h(k-1)) + imag(e(k)) );
        end

%         figure;
%         plot(1:N,real(h));
% 
%         figure;
%         plot(1:N,imag(h));

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

%         figure;
%         subplot(211);
%         semilogy(1:N, abs(h));
% 
%         subplot(212);
%         stem(1:N,abs(s-decision_matrix),"|");

        %% 7
%         error_matrix_magn = abs(error_matrix).^2;
        
%         for ind = 1:size( error_indices,1 )
%             if error_matrix_magn(error_indices(ind)) == 4
%                 err_bits = err_bits + 1;
%             elseif error_matrix_magn(error_indices(ind)) == 8
%                 err_bits = err_bits + 2;
%             end
%         end
        
        err_bits = err_bits + sum( abs( (error_matrix).^2 )/4 );

%         err_bits = err_bits + size( error_indices,1 );
    end
    snr_err_bits(snr_i) = err_bits;
%     disp(err_bits);
    snr_BER(snr_i) = err_bits/(N*packet_num);
%     disp(err_bits/(N*packet_num));
end

theoretical_BER = 1./(2*desired_snr_db);

figure;
plot(desired_snr_db, snr_BER);
hold on
plot(desired_snr_db, theoretical_BER);
xlabel('SNR','Interpreter','latex');
ylabel('BER','Interpreter','latex');
title('Experimental to theoretical BER');
legend('Experimental','Theoretical','Interpreter','latex');
legend show

figure;
semilogy(desired_snr_db, snr_BER);
hold on
semilogy(desired_snr_db, theoretical_BER);
xlabel('SNR','Interpreter','latex');
ylabel('BER','Interpreter','latex');
title('Experimental to theoretical BER');
legend('Experimental','Theoretical','Interpreter','latex');
legend show