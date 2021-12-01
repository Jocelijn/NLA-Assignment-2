diag_vals = [0:10];
A = diag(diag_vals);
v = 2.^diag_vals';%ones(10,1);
V = [v/norm(v)];
H = [];
[V,H] = ExtendArnoldi(A,V,H);
A*V- V*H
kAndResidual = []
for k = [1:9]
    sizeV = size(V)
    size_H = size(H)
    R = A*V(:,1:k)-V*H
    
    size_Hkk = size(H(1:k,1:k))
    [W,D] = eig(H(1:k,1:k))
    for j = [1:k]
        
        sizeW = size(W)
        eigenVec = V(1:k,1:k)*W(:,j);
        eigenVal = D(j,j);;
        
        r = norm(A*eigenVec );
        %newPoint = [k; r];
        %kAndResidual = [kAndResidual newPoint];
    end
    [V,H] = ExtendArnoldi(A,V,H);
end

plot(kAndResidual)

function [W,H2] = FilterAway(mu,V,H)
size_h = size(H);
k = size_h(2);
[Q,R] =  qr(H- mu*eye(k+1,k),0);
W = V*Q;

H2 = R*Q(1:k,1:k-1)+mu*eye(k,k-1);
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
    H
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

function [V,H] = ArnoldiMGS(A,v,k)
v = v/norm(v); 
V = v; H = [];
n = size(v);

for j=1:k
    w = A*v;   
    q = w;
    for l = 1:size(V,2)
        w = w-V(:,l)*(V(:,l)'*w);
    end
    g = norm(w);
    
    v = w/g;
    h = V'*q;
    
    V = [V v];
    H = [H,h; zeros(1,j-1),g];
end
end