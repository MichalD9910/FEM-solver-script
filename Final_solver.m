%"""Final_solver: Programm which solves heat conduction equations in two 
% dimension figure domain, using rectangular elements with 9 nodes"""
 
%__name__ = "Final_solver"
%__author__ = "Michal Dekarz"
%__version__ = "1.0"
%__status__ = "Complete" 


clc
close all
clear all

%___________Parameters___________%
%_Conduction in X direction_%
            k_x = 55; 

%_Conduction in Y direction_%
            k_y = 55;

%_Heat conduction factor_%
            h = 85;
            H = 10^8;

%_Heat flux q = dT/dx (Gęstość strumienia ciepła)_%
            q = 0; 
%_Convection_%
            T_inf = 23 + 273.15 ;

%_Temperature 1_%
            T1 = 30 + 271.15 ;

%_Temperature 2_%
            T2 = 100 + 271.15 ;


%___________Mesh___________%

%Defining parameters - size of cells in X and Y 
dx=1;
dy=1;
%Defining geometry points
g=[0 0;
    0 6; 
8 6;
8 3
4 3 
4 0; 
0 0];


%Amount of cells in X direction
podzialka1_x=floor((g(5,1)-g(1,1)))/dx;
podzialka2_x=floor((g(3,1)-g(1,1))/dx);

%Amount of cells in Y direction

podzialka1_y=floor((g(3,2)-g(1,2)))/dy;
podzialka2_y=floor((g(4,2)-g(1,2))/dy);



wezly=[];

for i=1:((g(2,2)-g(1,2))/dy+podzialka1_y+1)
    wezly_pionowo=unique([linspace(g(1,2),g(5,2),podzialka1_y+1) linspace(g(4,2),g(3,2),podzialka1_y+1)]);
    if i<=podzialka1_y
        wezly_poziomo=[linspace(g(1,1),g(5,1),podzialka2_x+1)];
    else
        wezly_poziomo=[linspace(g(1,1),g(4,1),2*podzialka2_x+1)];
    end
    wezly(size(wezly,1)+1:size(wezly,1)+length(wezly_poziomo),:)=[wezly_poziomo' wezly_pionowo(i)*ones(length(wezly_poziomo),1)];
    hold on;
end

%This figure shows the shape of domain
%figure(1)
%plot(wezly(:,1),wezly(:,2),'*')


%In this part points are getting connected to check the correctness of choice
tic
elementy = podzialka1_y*podzialka1_x+podzialka1_y*podzialka2_x
k=0;
l=0; 
m=2+podzialka2_x;
n=0;
u=0;
 
kk_global=zeros(size(wezly,1),size(wezly,1));
ka_global=zeros(size(wezly,1),size(wezly,1));
rb_global=zeros(size(wezly,1),1);


for j=1:elementy
    if j<=(podzialka1_y*podzialka1_x)
        if mod(j,rem(j,2))
            continue
           
        end
    
        if mod(j+k,(podzialka2_x)+1)==0
        k=k+podzialka2_x+2;
        end
  
        w1=j+k ; 
        w2=j+2+k;
        w3=j+4+2*podzialka2_x+k ; 
        w4=j+2+2*podzialka2_x+k ; 
        w5=j+1+k ;
        w6=j+k+podzialka2_x+3;
        w7=j+3+2*podzialka2_x+k;
        w8=j+k+podzialka2_x+1;
        w9=j+k+podzialka2_x+2;
        x1=wezly(w1,1); 
        y1=wezly(w1,2);
        x2=wezly(w2,1); 
        y2=wezly(w2,2);
        x3=wezly(w3,1); 
        y3=wezly(w3,2);
        x4=wezly(w4,1); 
        y4=wezly(w4,2);
        x5=wezly(w5,1); 
        y5=wezly(w5,2);
        x6=wezly(w6,1); 
        y6=wezly(w6,2);
        x7=wezly(w7,1); 
        y7=wezly(w7,2);
        x8=wezly(w8,1); 
        y8=wezly(w8,2);
        x9=wezly(w9,1); 
        y9=wezly(w9,2);        
     

    else 
     
    
        if mod(j,rem(j,2))
            continue
           
        end
        if n==podzialka2_x
     
            m=m+2*podzialka2_x+2;
            n=0;

        end
        w1=j+k+m+l; 
        w2=j+k+2+m+l;
        w3=j+4+4*podzialka2_x+k+m+l; 
        w4=j+2+4*podzialka2_x+k+m+l ;

        w5=j+1+k+m ;
        w6=j+k+m+2*podzialka2_x+3;
        w7=j+k+m+4*podzialka2_x+3;
        w8=j+k+m+2*podzialka2_x+1;
        w9=j+k+m+2*podzialka2_x+2;

        x1=wezly(w1,1); 
        y1=wezly(w1,2);
        x2=wezly(w2,1); 
        y2=wezly(w2,2);
        x3=wezly(w3,1); 
        y3=wezly(w3,2);
        x4=wezly(w4,1); 
        y4=wezly(w4,2);
        x5=wezly(w5,1); 
        y5=wezly(w5,2);
        x6=wezly(w6,1); 
        y6=wezly(w6,2);
        x7=wezly(w7,1); 
        y7=wezly(w7,2);
        x8=wezly(w8,1); 
        y8=wezly(w8,2);
        x9=wezly(w9,1); 
        y9=wezly(w9,2);
        n=n+1;
      


      
    end
    
    a1 = (x1 - x3)/2;
    b1 = (y1 - y4)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%In this part local function of shape is being solved .
syms u a s t b k_x k_y c alpha beta p q

A = [ones(9,1) [-a; 0; a; a; a; 0; -a; -a; 0] [-b; -b; -b; 0; b; b; b; 0; 0] [(-a)^2; (0)^2; (a)^2; (a)^2; (a)^2; (0)^2; (-a)^2; (-a)^2; (0)^2] [-a*(-b); 0*(-b); a*(-b); a*0; a*b; 0*b; -a*b; -a*0; 0*0] [(-b)^2; (-b)^2; (-b)^2; (0)^2; (b)^2; (b)^2; (b)^2; (0)^2; (0)^2] [((-a)^2)*(-b); ((0)^2)*(-b); ((a)^2)*(-b); ((a)^2)*0; ((a)^2)*b; ((0)^2)*b; ((-a)^2)*b; ((-a)^2)*0; ((0)^2)*0] [-a*(-b)^2; 0*(-b)^2; a*(-b)^2; a*(0)^2; a*(b)^2; 0*(b)^2; -a*(b)^2; -a*(0)^2; 0*(0)^2] [((-a)^2)*((-b)^2); ((0)^2)*((-b)^2); ((a)^2)*((-b)^2); ((a)^2)*((0)^2); ((a)^2)*((b)^2); ((0)^2)*((b)^2); ((-a)^2)*((b)^2); ((-a)^2)*((0)^2); ((0)^2)*((0)^2)]];

N = [1 s t s^2 s*t t^2 (s^2)*t s*(t^2) (s^2)*(t^2)] * inv(A); %f shape for element

B = [ diff(N,s,1); diff(N,t,1)];

C = [k_x 0; 0 k_y];

kk1 = B.'* C * B;
% kk=kk * B;

kk1 = int(kk1, t, -a, a);
kk1 = int(kk1, s, -b, b);

%kk1=subs(kk,a,0.25);
%kk1=subs(kk,b,0.25);
kk=subs(kk1,k_x,55);
kk=subs(kk,k_y,55); %macierz kk z podstawioną przewodnością

Nc12 = subs(N,s,c); 
Nc12 = subs(Nc12,t,-b); %f kształtu bok 1 - 2

Nc23 = subs(N,s,a);
Nc23 = subs(Nc23,t,c); %f kształtu bok 2 - 3

Nc34 = subs(N,s,-c);
Nc34 = subs(Nc34,t,b); %f kształtu bok 3 - 4

Nc41 = subs(N,s,-a);
Nc41 = subs(Nc41,t,-c); %f kształtu bok 4 - 1

k_p = -p * N.'* N;
k_p = int(k_p, t, -a, a);
k_p = int(k_p, s, -b, b);

r_q = q * N;
r_q = int(r_q, t, -a, a);
r_q = int(r_q, s, -b, b);

k_alpha12 = -alpha * Nc12.'* Nc12;
k_alpha12 = int(k_alpha12, c, -a, a);

k_alpha23 = -alpha * Nc23.'* Nc23;
k_alpha23 = int(k_alpha23, c, -b, b);

k_alpha34 = -alpha * Nc34.'* Nc34;
k_alpha34 = int(k_alpha34, c, -a, a);

k_alpha41 = -alpha * Nc41.'* Nc41;
k_alpha41 = int(k_alpha41, c, -b, b);

r_beta12 = beta * Nc12;
r_beta12 = int(r_beta12, c, -a, a);

r_beta23 = beta * Nc23;
r_beta23 = int(r_beta23, c, -b, b);

r_beta34 = beta * Nc34;
r_beta34 = int(r_beta34, c, -a, a);

r_beta41 = beta * Nc41;
r_beta41 = int(r_beta41, c, -b, b);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



     
    %plot([x1 x3],[y1 y3])
    %plot([x2 x4],[y2 y4])
    %plot([x5 x7],[y5 y7])
    %plot([x6 x8],[y6 y8])




%______________Local stiffness matrix______________%
      kk = subs(kk, a, a1);
      kk = subs(kk, b, b1);

      kk_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=kk_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+kk;


%Creating boundary conditions

%______________Boundary conditions______________%

      if y1==g(1,1) %boundary 100 degrees 
            alpha12= 85;
            beta12= -85 * (100+273.15);
            ka = subs(k_alpha12, a, a1);
            ka100 = subs(ka, alpha, alpha12);
            rb = subs(r_beta12, a, a1);
            rb = subs(rb, beta, beta12);
            rb100 = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka100;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb100;
      end
      if x2==g(3,2) %boundary 30 degrees 
            alpha23= 85;
            beta23= -85 * (30+273.15);
            ka = subs(k_alpha23, b, b1);
            ka30 = subs(ka, alpha, alpha23);
            rb = subs(r_beta23, b, b1);
            rb = subs(rb, beta, beta23);
            rb30 = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka30;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb30;
      end
      if x1==g(1,1) %boundary convection X
            alpha34= 85; 
            beta34= -85 * (23+273.15);
            ka = subs(k_alpha34, a, a1);
            ka_infX = subs(ka, alpha, alpha34);
            rb = subs(r_beta34, a, a1);
            rb = subs(rb, beta, beta34);
            rb_infX = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka_infX;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb_infX;
      end
      if y3==g(2,2) %boundary convection Y
            alpha41= 85;
            beta41= -85 * (23+273.15);
            ka = subs(k_alpha41, b, b1);
            ka_infY = subs(ka, alpha, alpha41);
            rb = subs(r_beta41, b, b1);
            rb = subs(rb, beta, beta41);
            rb_infY = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka_infY;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb_infY;
      end
      if y1==g(4,2) && x1>=g(5,1) %boundary flux Y
            alpha23= 0; 
            beta23= 0;
            ka = subs(k_alpha23, b, b1);
            ka_qY = subs(ka, alpha, alpha23);
            rb = subs(r_beta23, b, b1);
            rb = subs(rb, beta, beta23);
            rb_qY = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka_qY;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb_qY;

      end
      if x1==g(5,1) && y1<=g(4,2) %boundary flux X
            alpha12= 0; 
            beta12= 0;
            ka = subs(k_alpha12, a, a1);
            ka_qX = subs(ka, alpha, alpha12);
            rb = subs(r_beta12, a, a1);
            rb = subs(rb, beta, beta12);
            rb_qX = rb.';
            ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])=ka_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[w1 w5 w2 w6 w3 w7 w4 w8 w9])+ka_qX;
            rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])=rb_global([w1 w5 w2 w6 w3 w7 w4 w8 w9],[1])+rb_qX;

      end

      
end

%Visualization of the results

ui = inv(kk_global+ka_global)*rb_global;
ui=ui-273.15; 
figure(2)
[x,~,xi]=unique(wezly(:,1)); 
[y,~,yi]=unique(wezly(:,2));
subs=[xi yi];
M=accumarray(subs,ui,[],[],NaN);
[X,Y]=ndgrid(x,y);
surf(X,Y,M)
view(0,90),xlabel('x [cm]'),ylabel('y [cm]')
title("Heat graph in a flat model for the size of the grid:",num2str(dx)+" mm")
hc=colorbar;
colormap turbo; 
title(hc,'^\circ{C}');
saveas(gcf, [num2str(dx/2) '.jpg']);

%_________Temperature at boundary convection_________%

temp_x = M(1,:);
temp_y = M(2:end,end)';
temp_sr = mean(horzcat(temp_x,temp_y))
toc
