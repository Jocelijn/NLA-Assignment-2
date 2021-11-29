
%Todo: test
A = diag([0:10]);
v = ones(11,1);

LanczosMGS(A,v,4)

function T = abToT(a,b)
T =  diag(a,0) + diag(b,-1)+diag(b,1);
j = size(T,1);
z = zeros(j-1,1)';
T = [T; z b(end)];
end

function  LanczosMGS(A,v,k)
    v = v/norm(v); V = v;
    w = A*v; a(1) = v'*w;
    w = w - v*a(1); b(2) = norm(w);
    v = w/b(2); V = [V v];
    v_j = v;
    
    for j=2:k 
        Av_j = A*v_j;
        v_jmn1 = V(:,j-1);
        
        v_jp1 = Av_j;
        %orthogonalize v_jp1 with v_j-1:
        v_jp1 = v_jp1 - v_jmn1*(v_jmn1'*v_jp1);
        %orthogonalize v_jp1 with previous:
        v_jp1 = v_jp1 - v_j*(v_j'*v_jp1);
        

        %normalize:
        b(j+1) = norm(v_jp1);
        v_jp1 = v_jp1/norm(v_jp1);
        
        %we need that AV_j = B(j+1)V_j+1 + B(j)V_j-1 + A(j)V_j
        
        a(j) = Av_j'*v_j;
        v_j = v_jp1;
        V = [V v_jp1]
        
        T_j = abToT(a,b(2:j))
        
        [E,D] = eigs(T_j(1:j,1:j))
        u = V*E;
        residuals = A*u - D*u
        result = [vecnorm(residuals);diag(D)];
        result = [j;T_j(j+1,j) result]
        
    end
    size_a = size(a)
    size_b = size(b)
    T= abToT(a,b(2:k));
    
end