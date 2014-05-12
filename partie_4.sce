
// Question 13
function [u_next]= rk(u_n, t_n, f, dt)
    k1 = f(t_n, u_n)
    k2 = f(tn + dt * 0.5, un + dt * 0 .5 * k1)
    u_next = u_n + dt * 0.5 * (k1 + k2)
endfunction

//Question 14
functio
