
%%%---Input Parameter---%%%

E=eye(8);
e=exp(1);
rou=1.1;
Q2=zeros(8,1); 
x=zeros(8,1);             
x3=zeros(8,1);           
v2=zeros(8,1);
z2=zeros(8,1); 
FUN_it=zeros(1,100);
warning('off')


%%%---Signal Model---%%%
 H=[0.24,0.97,0.09,0.33,1.51,0.30,1.08,0.26;
    0.61,0.46,0.30,0.29,0.21,0.16,0.36,0.18;
    1.27,0.21,0.18,1.22,0.06,1.01,0.24,0.18;
    5.94,0.46,0.13,0.22,0.08,0.53,0.16,0.54;
    1.67,0.48,0.91,0.21,0.05,2.32,1.41,1.11;
    0.03,0.16,0.36,0.75,0.07,0.25,0.16,0.22;
    1.16,0.08,1.29,0.21,3.06,1.39,0.08,0.40;
    0.09,0.92,1.46,0.04,0.93,0.44,0.62,0.15];
 G=E;
 ZE=zeros(1,8);                      
 G(randperm(8,1),:)=ZE;
 H2=G*H;                      %channel 
 
 x=[0.25 0 0 0.5 0 0.75 0 0]';
 x=x/norm(x,2);     %transmitted signal
      
 y=H2*x;            
 y=awgn(y,15);      %received signal

         


  %%%---ADMM Algorithm---%%%
  L1=0.5;
  L2=0.5;
   
 v2=(H-H2)*x;  
 for k=1:150

     x3=(inv(2*(H2'*H2)+rou*E))*(rou*z2+2*(H2'*(y-v2))+Q2);        %x-parameters
     for t1=1:8
     if x3(t1,1)-(Q2(t1,1)/rou)>L1/rou             %z-parameters 
         z2(t1,1)=x3(t1,1)-(Q2(t1,1)/rou)-L1/rou;
     elseif x3(t1,1)-(Q2(t1,1)/rou)<-L1/rou
         z2(t1,1)=x3(t1,1)-(Q2(t1,1)/rou)+L1/rou;
     else z2(t1,1)=0;
     end 
     end
     for t2=1:8                                 %v-parameters
     if y(t2,1)-H2(t2,:)*x3>L2/2             
         v2(t2,1)=y(t2,1)-H2(t2,:)*x3-L2/2; 
     elseif y(t2,1)-H2(t2,:)*x3<-L2/2 
         v2(t2,1)=y(t2,1)-H2(t2,:)*x3+L2/2; 
     else v2(t2,1)=0;
     end
     end                                
      Q2=Q2+rou*(z2-x3);   %Q-parameters
     FUN_it(1,k)=norm(y-H2*x3-Q2,2)^2+0.5*norm(v2,1)+0.5*norm(z2,1);     %Objective Function
 end
  
  
  
  %%%---Output---%%%
  MSE=sum(x3-x).^2/8                      %MSE
