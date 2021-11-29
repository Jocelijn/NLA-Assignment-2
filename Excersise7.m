A = diag([1:10]);
v = ones(10,1);
k = 7;
[V,H] = ArnoldiMGS(A,v,k);
[W,H2] = FilterAway(2,V,H);

R = A*W(:,1:k-1)-W*H2
e = eig(H2(1:k-1,1:k-1))

[V2,H3] = ExtendArnoldi(A,V,H);
R = A*V2(:,1:k+1)-V2*H3
OrthoTest = V2*V2'
ListRitzData(H3)


function [W,H2] = FilterAway(mu,V,H)
size_h = size(H);
k = size_h(2)
[Q,R] =  qr(H- mu*eye(k+1,k),0)
size_q = size(Q)
size_r = size(R)
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
    H2 = [H,h; zeros(1,size(V,2)-1),g];
end

function L = ListRitzData(H)
    size_H = size(H)
    e = eig(H(1:size_H(2),1:size_H(2)));
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