function f = prior_pdf_twocomp(theta_trans)

	% function to compute prior density for mixture example

if (theta_trans(1) > theta_trans(2))
    f = 0;
else
    f = normpdf(theta_trans(1),0,1)*normpdf(theta_trans(2),0,1)*exppdf(exp(theta_trans(3)), 1)*exp(theta_trans(3))*exppdf(exp(theta_trans(4)), 1)*exp(theta_trans(4))*unifpdf(1/(1+exp(-theta_trans(5))),0,1)*exp(-theta_trans(5))/(1 + exp(-theta_trans(5)))^2;
end

end