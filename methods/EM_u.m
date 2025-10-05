function u = EM_u(A,sigma, u_ini,lanm_vec,r_vec, pho,b_up,b_low)

u = u_ini;
y_u =u;
inv_mat = inv(A'*A/sigma^2 + pho*eye(size(A,2)));
c = inv_mat*(lanm_vec+pho* r_vec);
t=1;
for i=1:60
    u_pre = u;

    a = (b_up - A*y_u)/sigma;
    b = (b_low - A*y_u)/sigma;
    pdf_diff = normpdf(a) -normpdf(b);
    cdf_diff = max(normcdf(a) -normcdf(b),1e-3);
    y = A* y_u - (pdf_diff./cdf_diff)*sigma;
    u = inv_mat*(A'*y)+c;
    
    t1=(1+sqrt(1+4*t^2))/2;
    y_u = u +(t-1)/t1*(u-u_pre);
    t=t1;
    % norm(u-u_pre,'fro')
    % if norm(u-u_pre,'fro')<1e-3
    %     break
    % end
end
