clear;
clc
n=50;%设置n
tic
syms b;
for i=1:2*n
syms (['x',num2str(i)]) ;
end
f=0;
for i=1:n
eval(['f','=','f','+','(',num2str(1),'-','x',num2str(2*i-1),')','^2','+','10*','(','x',num2str(2*i),'-','(','x',num2str(2*i-1),')','^2',')','^2',';']);
end
g=[];
for i=1:2*n
    t=diff(f,eval(['x',num2str(i)]));
    g=[g;t];
end
v=[];
for i=1:2*n
t=eval(['x',num2str(i)]);
v=[v,t];
end
G=hessian(f,v);
for i=1:2*n
eval(['x',num2str(i),'=','0',';']);
end
xian=cell(1);
xian{1,1}.dt=1;
xian{1,1}.x=zeros(2*n,1);
for i=1:2*n
    xian{1,1}.x(i,1)=double(subs(['x',num2str(i)]));
end
xian{1,1}.g=double(subs(g));
xian{1,1}.G=double(subs(G));
xian{1,1}.r=xian{1,1}.g;
xian{1,1}.f=0;
xian{1,1}.q=0;
xian{1,1}.rou=0;
k=1;
while(1)
   if (norm(xian{1,k}.g)<1e-5)
        disp(['信赖域法在迭代',num2str(k),'次终止']);
       break;
   end
   disp(['信赖域法第',num2str(k),'次迭代',':']);
  st=cell(1);
  st{1,1}.bt=0;
  st{1,1}.r=xian{1,k}.r;
  st{1,1}.p=-xian{1,k}.g;
  st{1,1}.x=zeros(2*n,1);
  for j=1:100
    if(norm(xian{1,k}.r)<=1e-5)
      xian{1,k}.s= zeros(2*n,1);
      disp(['子问题在迭代',num2str(j),'次满足停止测试']);
 break;
    end 
      if (st{1,j}.p'*xian{1,k}.G*st{1,j}.p<=0)
          t=double(solve((st{1,j}.x+b*st{1,j}.p).'*(st{1,j}.x+b*st{1,j}.p)==xian{1,k}.dt^2));
         s1=st{1,j}.x+t(1)*st{1,j}.p;
          s2=st{1,j}.x+t(2)*st{1,j}.p;
         if ((xian{1,k}.g'*s1+0.5*s1'*xian{1,k}.G*s1)>(xian{1,k}.g'*s2+0.5*s2'*xian{1,k}.G*s2))
           xian{1,k}.s= st{1,j}.x+t(2)*st{1,j}.p;
         else
            xian{1,k}.s= st{1,j}.x+t(1)*st{1,j}.p; 
         end
           disp(['子问题在迭代',num2str(j),'次遇到非正曲率']);
           break;
      end
      a=(st{1,j}.r'*st{1,j}.r)/(st{1,j}.p'*xian{1,k}.G*st{1,j}.p);
       st{1,j+1}.x=st{1,j}.x+a*st{1,j}.p;
      if(norm(st{1,j+1}.x)>=xian{1,k}.dt)
           t=double(solve((st{1,j}.x+b*st{1,j}.p).'*(st{1,j}.x+b*st{1,j}.p)==xian{1,k}.dt^2));
           s1=st{1,j}.x+t(1)*st{1,j}.p;
          s2=st{1,j}.x+t(2)*st{1,j}.p;
         if ((xian{1,k}.g'*s1+0.5*s1'*xian{1,k}.G*s1)>(xian{1,k}.g'*s2+0.5*s2'*xian{1,k}.G*s2))
           xian{1,k}.s= st{1,j}.x+t(2)*st{1,j}.p;
         else
            xian{1,k}.s= st{1,j}.x+t(1)*st{1,j}.p; 
         end
           disp(['子问题在迭代',num2str(j),'次达到信赖域边界']);
           break;
      end
      st{1,j+1}.r=st{1,j}.r+a*xian{1,k}.G*st{1,j}.p;
      if(norm(st{1,j+1}.r)<(1e-5*norm(xian{1,k}.r)))
           xian{1,k}.s=st{1,j+1}.x;
           disp(['子问题在迭代',num2str(j),'次满足停止测试']);
           break;
      end
      st{1,j+1}.bt=st{1,j+1}.r'*st{1,j+1}.r/(st{1,j}.r'*st{1,j}.r);
      st{1,j+1}.p=-st{1,j+1}.r+ st{1,j+1}.bt*st{1,j}.p;
  end
   if j==100
       xian{1,k}.s=zeros(20,1); 
       disp(['子问题在迭代',num2str(j),'次满足停止测试']);
   end
    xian{1,k+1}.x=xian{1,k}.x+xian{1,k}.s;
for i=1:2*n
eval(['x',num2str(i),'=',num2str(xian{1,k}.x(i,1)),';']);
end
f1=double(subs(f));    
for i=1:2*n
eval(['x',num2str(i),'=',num2str(xian{1,k+1}.x(i,1)),';']);
end
f2=double(subs(f));
xian{1,k}.f=f1-f2;
xian{1,k}.q=-xian{1,k}.g'*xian{1,k}.s-0.5*xian{1,k}.s'*xian{1,k}.G*xian{1,k}.s;
xian{1,k}.rou=xian{1,k}.f/xian{1,k}.q;
if (xian{1,k}.rou<0.25)
    xian{1,k+1}.dt=norm(xian{1,k}.s)/4;
elseif (xian{1,k}.rou>0.75&&norm(xian{1,k}.s)==xian{1,k}.dt)
    xian{1,k+1}.dt=2*xian{1,k}.dt;
else 
    xian{1,k+1}.dt=xian{1,k}.dt;
end
if(xian{1,k}.rou<=0)
    xian{1,k+1}.x=xian{1,k}.x;
else
    xian{1,k+1}.x=xian{1,k}.x+xian{1,k}.s;
end
for i=1:2*n
eval(['x',num2str(i),'=',num2str(xian{1,k+1}.x(i,1)),';']);
end
xian{1,k+1}.g=double(subs(g));
xian{1,k+1}.G=double(subs(G));
xian{1,k+1}.r=xian{1,k+1}.g;
if k==10
    break;
end
k=k+1;
end
toc





    

  
      
      
      




