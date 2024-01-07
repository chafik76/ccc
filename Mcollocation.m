function f=Mcollocation(n);
x=nodes(n);
%---------Number of the nodes of quadrature rule used in evaluating the integrals-------------------------------------------
m=100;
%---------------------------------The weight--------------------------------------------------------
lambda=[];
for i=1:n,
 lambda(i)=pi/n;
% lambda(i)=(pi/(n+1))*sin(i*pi/(n+1))^2;
end
%------------------------------------Function weight-----------------------------------------------------
w=@(x)(1-x).^(-(1/2)).*(1+x).^(-(1/2));
u=@(x)(1-x).^((0/3)).*(1+x).^((0/3));
%------------------------------------Matrices of the linear system-----------------------------------------------------
 A=zeros(n,n);
 S=zeros(n,n);
 D=zeros(n,n);
  for i=1:n,
      for j=1:n,
          if i==j,
            D(i,j)=u(x(i));
              else
              D(i,j)=0;
          end
      end
 end
%--------------------------------------------A-------------------------------------------------------------
 % for i=1:n,
      for j=1:n,
            A(:,j)=qua(@(z)kernel(z,x).*Lag(j,n,z),m); 
      % end
 end
%-----------------------------------------C--------------------------------------------------------------------
       for i=1:n,
              for j=1:n,
                        S(i,j)=qua(@(h)(kernel(h,x(i)).*qua(@(z)kernel(z,h).*Lag(j,n,z),m)),m);
             end
       end
%------------------------------------------The vectors of the linear system------------------------------------------------------
 b=zeros(n,1);
 d=zeros(n,1);
%------------------------------------------b------------------------------------------------------
  for i=1:n,
     b(i)=qua(@(s)kernel(s,x(i)).*smem(s),m);
  end
%------------------------------------------d------------------------------------------------------
  for i=1:n,
       d(i)=smem(x(i));
  end
%-------------------------Matrix and right hand side of linear system and solution----------------------------------------
  M1=D*(eye(n)-A)*inv(D);
  a1=M1\(D*d);  


R=D*(d+b-A*d);
  M=D*(eye(n)-S+A*A-A)*inv(D);
  a=M\R;

      y=0.1;

v=0;
% ----------------------------collocation---------------------------
    ss=-v;
  for j=1:n,
            ss=ss+(a1(j)/u(x(j))).*Lag(j,n,y);
  end
%  ----------------------------itertaed collocation---------------------------

  
  ss1=smem(y)-v;
   for j=1:n,
            ss1=ss1+(a1(j)/u(x(j)))*qua(@(z)kernel(z,y).*Lag(j,n,z),m);
  end











  %-------------------------Modified collocation----------------------------------------
    s=smem(y)-v;
  for j=1:n,
            s=s+(a(j)/u(x(j))-smem(x(j))).*Lag(j,n,y)+(a(j)/u(x(j)))*qua(@(z)kernel(z,y).*Lag(j,n,z),m);
  end
  s1=0;
 for i=1:n,
   for j=1:n,
            s1=s1-(a(j)/u(x(j)))*A(i,j)*Lag(i,n,y);
   end
  end
% -----------------------------------------Iterated modified collocation----------------------------------------
    s4=smem(y)+qua(@(s)kernel(s,y).*smem(s),m)-v;
  for j=1:n,
            s4=s4+((1/u(x(j)))*a(j)-smem(x(j))).*qua(@(s)kernel(s,y).*Lag(j,n,s),m)+...
              (a(j)/u(x(j)))*qua(@(h)(kernel(h,y).*qua(@(z)kernel(z,h).*Lag(j,n,z),m)),m);
  end
   s5=0;
  for i=1:n,
   for j=1:n,
            s5=s5-(1/u(x(j)))*a(j)*A(i,j)*qua(@(z)kernel(z,y).*Lag(i,n,z),m);
   end
  end
%    -------------------------Infinite norms of the four errors----------------------------------------
      f=[max(abs(ss.*u(y))) max(abs(ss1.*u(y))) max(abs((s+s1).*u(y)))  max(abs((s4+s5).*u(y)))];
    % f=max(abs((ss1).*u(y)));
    % f=cond(M);
%     f=max(abs((s4+s5).*u(y)));
% plot(y,s+s1);
%    hold on 
%     plot(y,s4+s5);
%     hold on 
%     plot (y,solex(y));
%  plot(y,s+s1-solex(y));
%     hold on 
%      plot(y,s4+s5-solex(y));
%     hold on 

%  f=cond(M);
 % % %------------------------------------------------------------------------------------------------------------------
% % % -----------------------------------------------------------------------------------------------------------------
% -------------------------------------------Zeros of jacobi polynomials-----------------------------------------------
function f=nodes(n);
x=[];
for i=1:n,
     x(i)=cos((2*i-1)*pi/(2*n));
%      x(i)=cos(i*pi/(n+1));
end
f=x;
% -------------------------------------------Lagrange polynomials-----------------------------------------------
 function f=Lag(i,n,y);
     x=nodes(n);
p=1;
for j=1:n,
    p=p.*((i==j)+(i~=j).*(y-x(j)))./((i==j)+(x(i)-x(j)));
end
f=p;
function K=kernel(x,y);
  % K=cos(y).*abs(x);
%  K=abs(x).^(5/2)./(y+2);
K=exp(y.*(x.^2+2))/30;
% function K=kernel2(x,y);
%  K=cos(y).*abs(x)*0.43816243616563694*pi;
%  K=exp(y.*(X.^2+2))/30;
% ------------------------------------Exacte solution--------------------------------------------------------------------------
% function f=solex(y);
%   f=abs(y).^(7/2);
%  f=abs(y).^(5/2);
% ------------------------------------second membre--------------------------------------------------------------------------
function f=smem(y);
  % f=solex(y)-(sqrt(pi)*cos(y)*gamma(11/4))./gamma(13/4);
  f=sign(y).*(1-y).^(1/3);
% -------------------------------------Averaged rule for computing the integrals----------------------------------------------------
function f=qua(g,m);
    s=(pi/(4*m))*(g(-1)+g(1))-g(-1)*(pi/(2*m));
    for i=1:m,
        s=s+(pi/(2*m))*(g(cos((2*i-1)*pi/(2*m)))+g(cos((m-i+1)*pi/(m))));
    end
    f=s;
    
    
   



function f=slx(n);
x=nodes(n);
%---------Number of the nodes of quadrature rule used in evaluating the integrals-------------------------------------------
m=n;
%---------------------------------The weight--------------------------------------------------------
lambda=[];
for i=1:n,
 lambda(i)=pi/n;
% lambda(i)=(pi/(n+1))*sin(i*pi/(n+1))^2;
end
%------------------------------------Function weight-----------------------------------------------------
w=@(x)(1-x).^(-(1/2)).*(1+x).^(-(1/2));
u=@(x)(1-x).^((0/3)).*(1+x).^((0/3));
%------------------------------------Matrices of the linear system-----------------------------------------------------
 A=zeros(n,n);
 S=zeros(n,n);
 D=zeros(n,n);
  for i=1:n,
      for j=1:n,
          if i==j,
            D(i,j)=u(x(i));
              else
              D(i,j)=0;
          end
      end
 end
%--------------------------------------------A-------------------------------------------------------------
 % for i=1:n,
      for j=1:n,
            A(:,j)=qua(@(z)kernel(z,x).*Lag(j,n,z),m); 
      % end
      j
 end
%-----------------------------------------C--------------------------------------------------------------------
 d=zeros(n,1);
%------------------------------------------d------------------------------------------------------
  for i=1:n,
       d(i)=smem(x(i));
  end
%-------------------------Matrix and right hand side of linear system and solution----------------------------------------
  M1=D*(eye(n)-A)*inv(D);
  a1=M1\(D*d);  


      y=-1:1/10:1;


%     ss=0;
%   for j=1:n,
%             ss=ss+(a1(j)/u(x(j))).*Lag(j,n,y);
%   end
% f=ss;
  ss1=smem(y);
   for j=1:n,
            ss1=ss1+(a1(j)/u(x(j)))*qua(@(z)kernel(z,y).*Lag(j,n,z),m);
   end
   f=ss1;
