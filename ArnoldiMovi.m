dimension = 50;
G = gallery('grcar',dimension);
v = ones(dimension,1);
% eigenV = eig(G)%LanczosCGS(G,v,50);
% size(eigenV)
% plot(eigenV,'*r')

h = figure;
eigenV_Ref = eig(G);
plot(eigenV_Ref,'*r');
axis tight manual
ax = gca;
ax.NextPlot = 'replacechildren';


loops = dimension;
M(loops) = struct('cdata',[],'colormap',[]);
h.Visible = 'off';
for j = 1:loops
    j
    
    eigenV = ArnoldiRitz(G,v,j);
    size(eigenV) 
    plot(eigenV_Ref,'*r');
    hold on;
    plot(eigenV,'*b');
    hold off;
    drawnow
    axis([-1,3, -4, 4])
    M(j) = getframe;
end

h.Visible = 'on';
movie(M);

function RitzValues = ArnoldiRitz(A,v,k)
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

H_kk = H(1:k,1:k);
RitzValues = eig(H_kk);
end