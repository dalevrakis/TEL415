clear all;
close all;

N = 2^6;
m = 100;
L = 3;
l = 0;
k = 2;
packets = 1000;
c = zeros(k,N);
h = zeros(k,L);
dec = zeros(1,m);
snrdb = [0:0.1:16];
tic
for SNR = 0:2:16
    errors1 = 0;
    for p= 1:packets
        s = ones([k m]);
        for i = 1:k
            c(i,:) = 1/sqrt(N)*sign(randn(1,N));
            h(i,:) = sqrt(1/(L*2))*(randn(1,L)+1i*randn(1,L));
            for j=1:m
                random = randi(2);
                if random==1
                    s(i,j)= 1;
                elseif random==2
                    s(i,j)= -1;
                end
            end

        end

        sigma = 1/(10^(SNR/10));
        % y = zeros(N+L-1,1);
        % for i = 1:k
        %     w = sqrt(sigma/2)*sign(randn(1,N));
        %     y(i,:) = c(i,:).*s;
        % end 
        for i = 1:m
            y = zeros(N+L-1,1);
            
            for j = 1:k
                for l = 0:L-1
                    y = y + h(j,l+1)*s(j,i)*a_l(c(j,:),l,L);
                end
            end

            y = y + sqrt(sigma/2)*(randn(N+L-1,1)+1i*randn(N+L-1,1));

            norm_h = norm(h(1,:),2);
            hs = zeros(k,L);
            for j = 1:k
                for l = 0:L-1
                   hs(j,l+1) = sum(conj(h(j,l+1))/norm_h.*y.*a_l(c(j,:),l,L));
                end
            end
            r = sum(hs(1,:));

            dec(i) = sign(real(r));
        end
        errors1 = errors1 + sum(abs(s(1,:)-dec).^2)/4;

    end
    BER1(SNR/2+1)= errors1/(m*packets)
end
toc
figure()
semilogy(0:2:SNR,BER1)
hold on
semilogy(snrdb, 1./snrdb.^1)
semilogy(snrdb, 1./snrdb.^2)
semilogy(snrdb, 1./snrdb.^3)
xlabel("SNR_{db}")
ylabel("BER")
% legend('MRC','TB','fontsize',15,'interpreter','latex','Location','southwest')
% exportgraphics(gcf,"asurmates27.png",'Resolution',300);

function a = a_l(c,l,L)
    a = [zeros(1,l) c zeros(1,L-1-l)]';
end