function VE=VEFI_to_VE(VEFI)

phi=VEFI(1);
theta=VEFI(2);
psi=VEFI(3);
phip=VEFI(4);
thetap=VEFI(5);
psip=VEFI(6);

VE(4:7) = QuaternionfromEuler(phi,theta,psi); 
VE(1:3)=[phip-psip*sin(theta)  ;
         thetap*cos(phi)+psip*cos(theta)*sin(phi);
            -thetap*sin(phi)+psip*cos(theta)*cos(phi)  ];% Ya que es solo una operacion aislada utilizaremos la ecuaciones con rel trigonometricas


end