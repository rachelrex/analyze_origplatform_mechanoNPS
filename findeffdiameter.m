function [Deffvalue] = findeffdiameter(d,L,I)

%d= known bead diameter
%L=length of channel 
%I=deltaI/I 
syms Deff
left=(d^3)/((Deff^2)*L);
right= 1/(1-(0.8*((d/Deff)^3)));
eqn=(left*right)-I;
Deffvalue= vpasolve(eqn,Deff);
