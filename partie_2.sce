//question 5 :
exec('partie_1.sce')
//avec u=(e,r)
function coord(M,i,j,u)
    M(:,:,i,j) = u
endfunction
    
function [matrice]=gen_grille(n,t,u)
    h=1/(n-1)
    matrice=hypermat([2 1 n n])
    for i=1:n
        for j=1:n
            coord(matrice,i,j,u)
        end
    end
    return matrice
endfunction

//question 6
// Permet d'ecrire un bloc matrice dans la matrice L
function M=write_T(M,T,x,y,n)
    for i=1:1:n
        for j=1:1:n
            M(x-1+i,y-1+j)= T(i,j)
        end
    end
    return M
endfunction

function L=laplacien(D,n)
    h=1/(n-1)
    T_base = zeros(n,n)
    for i=2:1:(n-1)
            T_base(i,i-1)=1
            T_base(i,i+1)=1
    end
    T_base(1,2)=1
    T_base(n,n-1) = 1
    
    T0 = T_base
    for i=2:1:(n-1)
        T0(i,i)=-3
    end
    T0(1,1)=-2
    T0(n,n)=-2
    T0=1/(h*h)*T0
    TN=T0
        
    TK = T_base
    for i=2:1:(n-1)
        TK(i,i)=-4
    end
    TK(1,1)=-3
    TK(n,n)=-3
    TK=1/(h*h)*TK
    
    L = zeros(n*n,n*n)
    L=diag(1/h**2*ones(n**2-n,1),n) + diag(1/h**2*ones(n**2-n,1),-n)
    
    //ecriture des T et I
    for i=1:n:n*n
       L=write_T(L,TK,i,i,n)
    end
    L=write_T(L,T0,1,1,n)
    L=write_T(L,T0,n**2-n+1,n**2-n+1,n)
    L=D*L
    
endfunction

//question 7
function [Ls]=slaplacien(D,n)
    Ls=sparse(laplacien(D,n))
endfunction

//question8

function test_ls(n)
    B=ones(n*n,1)
    ls=slaplacien(10,n)
    ls*B
endfunction

function test_l(n)
    B=ones(n*n,1)
    l=laplacien(10,n)
    l*B
endfunction

//question 9
function y=F_2d(E,R)
    y=-k*E*(E-alpha1)*(E-1)-E*R
endfunction
function y=G_2d(E,R)
    y=(epsilon+mu1*R/(E+mu2))*(-R-k*E**2-alpha2-1)
endfunction

function y=L(E,A)
    y=A.*E'
endfunction

function [u]=FG(dt,un)
    u(1)=F_2d(un(1),un(2)) + Le_k
    u(2)=G_2d(un(1),un(2))
endfunction

function [E,R]=question9(n,t_max,E,R,H)
    k=8.0
    alpha1=0.1
    alpha2=0.1
    epsilon=0.01
    mu1=0.1
    mu2=0.1
    dt=1/(t_max)
    A = slaplacien(1,n)
    D=1
    
    
    e0=E(1,1,1)
    r0=R(1,1,1)
    vect_n=[e0;r0]
    
    x=[0:h:1]
    vect_E=[]
    Le=[]
    
    for i=1:1:45
        //calcule de L
        vect_E=[]
        for r=1:1:n
            vect_E=[vect_E,E(r,:,i)]
        end
        Le=A*vect_E'
        
        for j=1:1:n
            for k=1:1:n
                Le_k= Le(j*k)
                vect_n=euler1(vect_n,10,dt,H)
                E(j,k,i+1)=vect_n(1)
                R(j,k,i+1)=vect_n(2)
            end
        end
    end
    
endfunction



//question 10

function [u]=new_FG(dt,un)
    u(1)= un(1)*(2*%pi**2*D-1) + D*2*%pi**2*un(1)
    u(2)=-un(2)
endfunction

function y=u_a(u,v,t)
    y=cos(u*%pi)*cos(v*%pi)*exp(-t)
endfunction

function y=v_a(u,v,t)
    y=cos(2*u*%pi)*cos(2*v*%pi)*exp(-t)
endfunction

function [E,R]=init_ER(n)
    
    for i=1:n
        for j=1:n
            E(i,j,1)=u_a(i*h,j*h,0)
            R(i,j,1)=v_a(i*h,j*h,0)
        end
    end

endfunction

function plot_ER_3d(E,t)
    x=[0:h:1]
    for i=1:n
        for j=1:n
            z(i,j)=E(i,j,t)
        end
    end
    plot3d(x,x,z)
endfunction

function exemple(n,t_max)
    h=1/(n)
    [E,R]=init_ER(n)
    
   
    [E,R]=question9(n,t_max,E,R,new_FG)
    
    
    //plot_ER_3d(E,2)
    x=[0:h:1]
    y=[0:h:1]
    
    for i=1:n
        for j=1:n
            z(i,j)=u_a(i*h,j*h,2)
        end
    end
    //plot2d(x,E(1,:,2))
    //plot3d1(x,y,z)
    disp(E(:,:,1))
    disp(E(:,:,2))
    disp(E(:,:,3))
endfunction

    


function verif9()
    
endfunction

//question 11 :
function simulation2d(t_max,n)
     E=hypermat([n n t_max])
     R=hypermat([n n t_max])
     
     //initialisation de E et R 
        h=1/(n-1)
        for i=1:1:n
            if i*h < 0.5 then
                E(i,:,1)=2
                R(:,i,1)=0
            else
                E(i,:,1)=3
                R(:,i,1)=1
            end
        end
    
      

     [E,R]=question9(n,t_max,E,R,FG)
    plot_ER_3d(E,2)
    plot_ER_3d(R,1)
    
endfunction
//simulation2d(400,10)
exemple(3,2000)
