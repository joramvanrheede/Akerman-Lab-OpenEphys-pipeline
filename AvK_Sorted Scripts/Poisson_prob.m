function Out = Poisson_prob(x,lamda);
Out = ((lamda^x)*(exp(1)^(0-lamda)))/factorial(x);
end