function H=nak_m(m,Nr,Nt)
n=zeros(Nr,Nt);
for i=1:2*m
n=n+randn(Nr,Nt).^2;
end
n=n/(2*m);
phi=2*pi*rand(Nr,Nt);
H=(n.^0.5).*cos(phi)+1i*(n.^0.5).*sin(phi);