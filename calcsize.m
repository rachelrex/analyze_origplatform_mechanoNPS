function size=calcsize(deff,L,deltRR)
top=deltRR*(deff^3)*L;
bottom=deff+(0.8*L*deltRR);
size=(top/bottom)^(1/3);
%test