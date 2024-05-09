function [current_pt, stats] = ORN_ls_simple(F_func, J_func, sigma_0, mu, alpha, beta, init_pt, N_iter, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second-order optimsitc method with line search
% for solving the nonlinear equation F(z)=0
% with a monotone operator F
% Inputs:
%   F_func: the operator F(z)
%   J_func: the Jacobian DF(z)
%   sigma_0: the initial trial stepsize
%   mu: the strongly-convex parameter
%   alpha: line search parameter
%   beta: line search parameter
%   init_pt: initial point
%   N_iter: max. number of iterations
%   epsilon: error tolerance
%
% Outputs:
%   list_iter: the generated iterates {z_k}
%   list_dis: the displacement ||z_{k+1}-z_k||
%   list_eta: the stepsize
%   list_steps: the number of calls to the optimistic subsolver
%   flag: true if we stop before reaching the max. number of iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ruichen Jiang (rjiang@utexas.edu), UT Austin, 2024
%
% Modified by Joao M Pereira for compatibility with npe and hipnex
% experiments

% list_iter = zeros(length(init_pt),N_iter+1); % record the iterates
fcn_timer = tic;
list_dis = zeros(N_iter,1); % record the displacements
list_eta = zeros(N_iter,1); % record the stepsizes
list_steps = zeros(N_iter,1); % record the number of LS steps
it_times = zeros(N_iter, 1);
% list_iter(:,1) = init_pt;


F_norms = zeros(N_iter,1);


current_pt = init_pt;
v = 0*init_pt; % scaled residual v
% eta_up = sigma_0;
sigma = sigma_0;
F = F_func(current_pt);
n = length(F);

F_evals = 1;
DF_evals = 0;
linsolve_count = 0;
lap = toc(fcn_timer);
it_times(1) = lap;
F_norms(1) = norm(F);

for it = 1:N_iter
    % warmstart: the initial trial stepsize is the inadmissible one 
    % in the last iteration
    % sigma = eta_up;
    eta = sigma;
    J = J_func(current_pt);
    DF_evals = DF_evals + 1;
    flag_stepsize = false;
    step_counter = 0;

    while ~flag_stepsize
        if step_counter >0
            eta = beta*eta; % backtrack
        end
        
        direction = (eta*J+eye(n))\(-eta*F-v);
        linsolve_count = linsolve_count + 1;
        pt_new = current_pt+direction;
        approx_first = F+J*(pt_new-current_pt);
        dis = norm(direction);
        F_pt_new = F_func(pt_new);
        F_evals = F_evals + 1;
        res_new = F_pt_new-approx_first;
        flag_stepsize = eta*norm(res_new)<=1/2*alpha*dis;

        step_counter = step_counter+1;
    end
%         % record the admissible stepsize
%         eta_lo = eta;
%         pt_lo = pt_new;
%         F_lo = F_pt_new;
%         v_lo = res_new/(1/eta+mu);
%         dis_lo = dis;
    % list_iter(:,1+i_iter) = pt_new;
    list_dis(it) = dis;
    list_eta(it) = eta;
    list_steps(it) = step_counter;

    current_pt = pt_new; % update the iterate
    F = F_pt_new; % update the operator (gradient)
    v = res_new/(1/eta+mu); % update the scaled residual
    sigma = eta * sqrt(1+2*eta*mu) / beta;
    norm_F = norm(F_pt_new);
    F_norms(it+1) = norm_F;
    new_lap = toc(fcn_timer);
    it_times(it+1) = new_lap - lap;
    lap = new_lap;
            
    if norm(F_pt_new)<= epsilon || dis<=epsilon
        break
    end
end
% list_iter = list_iter(:,1:1+i_iter);
stats.list_dis  = list_dis(1:it);
stats.list_eta = list_eta(1:it);
stats.list_steps = list_steps(1:it);
stats.it_times = it_times(1:it+1);
stats.F_norms = F_norms(1:it+1);
stats.lambdas = list_eta(1:it);
stats.linsolve_count = linsolve_count;
stats.n_iters = it;
stats.F_evals = F_evals;
stats.DF_evals = DF_evals;

end
