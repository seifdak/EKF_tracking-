function df = fdjac(x,g,varargin)

fx = feval(g, x, varargin{:});
[F,C] = size(fx);
fh = zeros(F,C);
escolumna = 1;
if (C~=1)   
    escolumna=0;    
end

Lfx= length(fx);
Lx = length(x);
df = zeros(Lfx, Lx);
for j=1:Lx
    temp = x(j);
    h = 1e-8;
    x(j) = temp+h;
    fh = feval(g,x,varargin{:});
    x(j)=temp;
    if escolumna    
       df(:,j)=(fh-fx)/h;
    end
end


