function p = TauExact(Bounds,sigma)

% Specifying parameters
nx=80;                           
ny=80;                                  
niter=1000;                       
dx=2/(nx-1);                     
dy=2/(ny-1);                    
x=Bounds(1,1):dx:Bounds(1,2);                      
y=Bounds(2,1):dy:Bounds(2,2);                  
b=zeros(nx,ny);                        

% Initial Conditions
p=zeros(nx,ny);                 

%Boundary conditions
p(:,1)=0;
p(:,ny)=0;
p(1,:)=0;                  
p(nx,:)=0;

%Source term
b(:,:) = -2/(sigma^2);

i=2:nx-1;
j=2:ny-1;

for it=1:niter
    pn=p;
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(b(i,j)*dx^2*dy*2))/(2*(dx^2+dy^2));
    p(:,1)=0;
    p(:,ny)=0;
    p(1,:)=0;                  
    p(nx,:)=0;
end

h=surf(x,y,p','EdgeColor','none');       
xlabel('x')
ylabel('y')
zlabel('\tau')