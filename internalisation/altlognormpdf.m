function density = altlognormpdf(x,mu,sigma)

    d = makedist('Lognormal',mu,sigma);
    
    density = pdf(d,x + mean(d) - 1);

end