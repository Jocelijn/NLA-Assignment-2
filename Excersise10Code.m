diag_vals = [0+1i, 0-1i,1+2*1i, 1 - 2*1i, 3 - 3*1i, 3 + 3*1i,5,6,9,10];
A = diag(diag_vals);
v = ones(10,1);
V = [v/norm(v)];
H = [];
[V,H] = ExtendArnoldi(A,V,H);
[V,H] = ExtendArnoldi(A,V,H);
[V,H] = ExtendArnoldi(A,V,H);
[V,H] = ExtendArnoldi(A,V,H);
[V,H] = ExtendArnoldi(A,V,H);
[V,H] = ExtendArnoldi(A,V,H);
[W,H2] = FilterAway2(0+1i,V,H,v,A)

while true
   
    c = input("What do you want to do? chose ex, filter, or exit ","s")
   
    if strcmp(c,"filter")
        disp("The Ritz pairs and their resiudals are:")
        mu = ListRitzData(H);
        j = input("input the index of the eigenvalue you want to remove");
        
        [V,H] = FilterAway(mu(j),V,H);
    end
    
    if strcmp(c,"ex")
        [V2,H2] = ExtendArnoldi(A,V,H);
       V = V2;
       H = H2
       mu = ListRitzData(H);
       smallest_eigen = min(mu)
       
    end
    
    if strcmp(c,"exit")
        ListRitzData(H)
        return
    end
end


function [W,H2] = FilterAway(mu,V,H)
size_h = size(H);
k = size_h(2);
[Q,R] =  qr(H- mu*eye(k+1,k),0);
W = V*Q;

H2 = R*Q(1:k,1:k-1)+mu*eye(k,k-1);
end

function [W,H2] = FilterAway2(mu,V,H, oldStart,A)
size_h = size(H);
size_a = size(A)
n = size_a(2);
k = size_h(2)



V_k = V(:,1:k);
H_k_kmin1 = H(1:k,1:k-1);
H_kmin1_kmin2 = H(1:k-1,1:k-2);
V_kmin1 = V(:,1:k-1);
V_kmin2 = V(:,1:k-2);



P = V_k*H_k_kmin1*H_kmin1_kmin2 ;
P = P - (mu + conj(mu))*V_kmin1*H_kmin1_kmin2 ;
P = P + abs(mu)^2*V_kmin2;

size(P)
size([1;zeros(k-3,1)]')
newStart = P* ([1;zeros(k-3,1)])
newStart2 =(A- mu*eye(n))*(A- conj(mu)*eye(n))*oldStart %H*H - (mu + conj(mu))+ abs(mu)^2*eye(k+1,k);

[Q,R] =  qr(P,0);

size(Q)
size(V)
size(V(:,1:k-2))
W = (V(:,1:k-2)')*(Q);
size(W)
H2 = Q'*H*Q;%R*Q(1:k,1:k-1)+mu*eye(k,k-1);
end

function [V2,H2] = ExtendArnoldi(A,V,H)
    v = V(:,end);
    w = A*v;   
    q = w;
    for l = 1:size(V,2)
        w = w-V(:,l)*(V(:,l)'*w);
    end
    g = norm(w);
    
    v = w/g;
    h = V'*q;
    
    V2 = [V v];
    H2 = [H,h; zeros(1,size(V,2)-1),g];
end

function L = ListRitzData(H)
    size_H = size(H);
    [V,D] = eig(H(1:size_H(2),1:size_H(2)));
    k = size_H(2);
    h = H(k+1,k);
    ResultTable = [k; abs(h)];
    for j = 1 : k
        eigenVec = V(:,j);
        r = abs(h*eigenVec(k));
        eigenVal = D(j,j);
        
        newColumn = [eigenVal;r];
        ResultTable = [ResultTable newColumn];
    end
    display(ResultTable)
    L =  diag(D);
end