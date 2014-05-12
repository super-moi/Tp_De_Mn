//*************************************************************//
// question 1 :
//calcule d'un pas du probleme d'Euler
function [vect_suivant] = euler1(vect_n, t, dt, fonct)
        vect_suivant = vect_n + dt*fonct(t,vect_n)
endfunction
//*************************************************************//


//*************************************************************//
// question 2 :
//fonction exemple
function  y = g(t,u)
       y = -0.1 * t * u
endfunction

//fonction vérifiant la fonction euler1
function verif_euler1(f,u0)
    dt = 0.01
    n = 0.0
    u_n = u0
    y=[0.0:dt:5.0]
    t_n = 0.0
    index = 1
    while n < 5.0
        u_n = euler1(u_n, t_n, dt, f)
        y(index)=u_n
        n = n + dt
        t_n = t_n + dt
        index = index + 1
    end
    x=(0.0:dt:5.0)
    plot2d(x,[y',exp(-0.05*(x')**2)])
endfunction
//*************************************************************//

//*************************************************************//
//question 3 :
// fonction F et G du probleme 
function y=F(e,r)
      y=-k*e*(e-alpha1)*(e-1.0)-e*r
endfunction
function y=G(e,r)
    y=(epsilon+mu1*r/(e+mu2))*(-r-k*e*(e-alpha2-1))
endfunction

//systeme du schema d'Euler
function [K]=H(dt,u)
    K(1)=F(u(1),u(2))
    K(2)=G(u(1),u(2))
endfunction
   
function [sol]=question3()
    //initialisation des parametres
    n=50
    k=8.0
    alpha1=0.1
    alpha2=0.1
    epsilon=0.01
    mu1=0.1
    mu2=0.1
    e0=1.0
    r0=0.0
        
    dt=0.01
    x=[0]
    
    vect_n = [e0;r0]
    sol = vect_n
    
    index=1
    //résolution :
    for t=0:dt:n
        vect_n=euler1(vect_n,t,dt,H)
        sol=[sol,vect_n]
        x=[x,t] 
        dt=0.01
        index=index+1
    end
     
    plot2d(x,sol')
    return sol
    disp(sol)

endfunction
//*************************************************************//
