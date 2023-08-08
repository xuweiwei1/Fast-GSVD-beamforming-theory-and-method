function [U,V,X,psi,phi]=fgsvd(A,B)
[m,n]=size(A);
[p,~]=size(B);
H=[A;B];
[L,R]=qr(H,0);
L1=L(1:m,:);
L2=L(m+1:end,:);
q=size(L1,2);
if n==q
    if m<p
        [S,psi]=eig(L1'*L1);
        phi=diag(sqrt(ones(q,1)-diag(psi)));
        psi=diag(sqrt(diag(psi)));
        if q>m
        psi=psi(q-m+1:q,:);
        end
        if q>p
        phi=phi(1:p,:);
        end
    else
        [S,phi]=eig(L2'*L2);
        psi=diag(sqrt(ones(q,1)-diag(phi)));
        phi=diag(sqrt(diag(phi)));
        if q>p
        phi=phi(q-p+1:end,:);
        end
        if m<q
        psi=psi(1:m,:);
        end
    end
    X=S'*R; 
    U=L1*S/psi;
    V=L2*S/phi; 
else
        [U,~]=eig(L1*L1');
        [V,~]=eig(L2*L2');
        psi=[eye(m) zeros(m,q-m)];
        phi=[zeros(p,q-p) eye(p)];
        X=[U'*A;V'*B];
end
end
