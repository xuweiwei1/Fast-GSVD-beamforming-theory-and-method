function [s] = QPSK_mapper(q,Tc,S,K)
s=zeros(K*S,Tc/2);
s_real=0;
s_imag=0;
for i1=1:K*S,
    for i2=1:Tc/2,
        switch q(i1,2*i2-1),
            case 0,
                s_real=(sqrt(2))/2;
            case 1,
                s_real=-(sqrt(2))/2;
        end
        switch q(i1,2*i2),
            case 0,
                s_imag=(sqrt(2))/2;
            case 1,
                s_imag=-(sqrt(2))/2;
        end
       s(i1,i2)=(1/sqrt(S))*complex(s_real,s_imag);
    end
end