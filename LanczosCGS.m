function LanczosCGS(A,v,k) 
T = zeros(k+1,k);

%first iteration
v = v/norm(v);      V = v;
w = A*v;            a(1) = v'*w;
w = w - v*a(1);     b(2) = norm(w);
v = w/b(2);         V = [V v];

%forming of first column of T
t = [a(1) b(2) zeros(1,k-1)];
T(:,1) = t';

%solving eigenvalue problem
[y,mu] = eig(T(1,1));
Mu=[1];
Mu = [Mu, mu(1,1)];
r = [abs(T(2,1))*abs(y(1,1))];

%creating table (output)
Table=[Mu; T(k+1,k) r];
disp(Table)

%iterations 2 - k
for j=2:k
r = [];
%Lanczos iteration
w = A*v;
a(j) = v'*w;
w = w - b(j)*V(:,j-1)-a(j)*v;
b(j+1) = norm(w);
v = w/b(j+1);
V = [V v];  

%forming of rest of the columns of T
t = [zeros(1,j-2) b(j) a(j) b(j+1) zeros(1,k-j)];
T(:,j) = t';

%solving eigenvalue problem
[y,mu] = eig(T(1:j,1:j));
Mu=[j];
for i=1:j
    Mu = [Mu, mu(i,i)];
    r = [r,abs(T(i+1,i))*abs(y(i,i))];
end

%creating table (output)
Table=[Mu; T(k+1,k) r];
disp(Table)

end

for n=1:j
    Plot(j,n)=Table(1,n+1);
end

end

color = rand(k,3);
for t=1:k
    plot(t,Plot(t,1:t), '-*','MarkerEdgeColor', color(t,:));
    hold on;
end
