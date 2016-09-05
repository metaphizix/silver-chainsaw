clc
clear all
disp('Variables')
%Variables
x_0 = [1,0]' %
Delta_0 = 1
Delta_k = 1
n = 1/4

xk = x_0
%Hessian Matrix
syms x1 x2  
f = 10*(x2 - x1.^2).^2 + (x1-1).^2;
B_k = hessian(f)

%Hessian Matrix at x_0

B_k_Atx_0 = subs(B_k,{x1,x2},{x_0(1),x_0(2)});
disp('Hessian at x_0')
B_k_Atx_0 

gradf = gradient(f, [x1, x2,])


% Checking if Matrix is Positive Definite
    eig_A = eig(B_k_Atx_0);
    flag = 0;
    for i = 1:rank(B_k_Atx_0)
        if eig_A(i) <= 0 
        flag = 1;
        end
    end

    if flag == 1
        disp('the matrix is not positive definite')
    else
        disp('the matrix is positive definite')
    end

    
grad_f_xk = subs(gradf,{x1,x2},{xk(1),xk(2)})
norm(grad_f_xk)
disp('Trust Region Algorithm')
iterations = 0;

%while subs(gradf,{x1,x2},{xk(1),xk(2)}) ~= 0
%while norm(subs(gradf,{x1,x2},{xk(1),xk(2)})) <= 0.01    
%while true
while iterations < 100
    %------------Second Step of The Trust Region--------------------------
   
    B_k_Atxk = subs(B_k,{x1,x2},{xk(1),xk(2)});

    % Checking if Matrix is Positive Definite
    eig_A = eig(B_k_Atxk);
    flag = 0;
    for i = 1:rank(B_k_Atxk)
        if eig_A(i) <= 0 
        flag = 1;
        end
    end

    if flag == 1
        %disp('the matrix is not positive definite')
    else
        %disp('the matrix is positive definite')
    end

    

    %grad at x

    grad_f_xk = subs(gradf,{x1,x2},{xk(1),xk(2)})
    norm(grad_f_xk)
  
    
    % Calculating Tau
    if grad_f_xk'*B_k_Atxk*grad_f_xk  <= 0
        tau = 1;
    else
        tau = ((norm(grad_f_xk))^3)/(Delta_k*grad_f_xk'*B_k_Atxk*grad_f_xk);
    end

    % Calculating Steepest Descent Direction
    d_k_b = -(inv(B_k_Atxk))*grad_f_xk;
    
    %Verify if Norm of d_k_b is less than Delta
    %disp('Norm of d_k_b')
    %norm(d_k_b)
    %disp('Delta')
    Delta_k;


    d_k_u = -((grad_f_xk'*grad_f_xk)/(grad_f_xk'*B_k_Atxk*grad_f_xk))*grad_f_xk;

    if tau <= 1 && tau >=0
        d_k_tau = tau*d_k_u;
    elseif tau <= 2 && tau >=1
        d_k_tau = d_k_u + ((tau-1)*(d_k_b-d_k_u));
    end    
  
    % Calculating Actual Reduction and Prediction Reduction

    %Verify if Norm of d_k_tau is less than Delta
    %disp('Norm of d_k_tau')
    %norm(d_k_tau)
    %disp('Delta')
    Delta_k;

    %Quadratic Model

    f_xk = subs(f,{x1,x2},{xk(1),xk(2)});
    m_k_d = f_xk + (grad_f_xk'*d_k_tau)+(0.5*(d_k_tau'*B_k_Atxk*d_k_tau));
    xk_d = xk + d_k_tau;
    f_xk_d = subs(f,{x1,x2},{xk_d(1),xk_d(2)});
    p_k = (f_xk - f_xk_d)/(f_xk - m_k_d);
    
    if f_xk == m_k_d ,break ,end
       

    %------------Third Step of The Trust Region--------------------------
    if p_k > 0.25 
        Delta_k = 0.25*norm(d_k_tau);
    elseif p_k > 0.75 && norm(d_k_tau)== Delta_k    
        Delta_k = min(2*Delta_k,Delta);
    else
        Delta_k = Delta_k;
    end

    %------------Fourth Step of The Trust Region--------------------------
    if p_k > n 
        xk = xk_d;
        %continue
    else
        xk = xk;
    end
    %disp('x after iteration')
    
    iterations = iterations +1;
    
    xk ;
    %After This point the Algorith Diverges
    if norm(grad_f_xk)<=0.1172
        break
    end    
    
end    
disp('------------------------------------------------------------------')

disp('Norm of gradf')
StoppinCondition = norm(subs(gradf,{x1,x2},{xk(1),xk(2)}))
disp('Minimizer is')
xk

str = fprintf('Number of iterations = %d', iterations );