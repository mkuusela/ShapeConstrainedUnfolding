% Objective function for fmincon
function [fun,grad] = funCon(nuTilde,f)

fun = f'*nuTilde;
grad = f;