(* Created with the Wolfram Language : www.wolfram.com *)
{
 {(I*((4*M*(a0 + 2*MBH)*(-2*MBH + r)^2)/
       (E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 
           2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*r*(a0 + r)^3) + 
      r^2*\[Omega]^2)*H0[r])/
    (r^2*(\[Omega] - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2)*\[Omega])/r)) + 
   (2*(a0^3*MBH + 3*a0^2*MBH*r + MBH*r^3 + M*r*(8*MBH^2 - 6*MBH*r + r^2) + 
      a0*MBH*(4*M*MBH - 2*M*r + 3*r^2))*H1[r])/(r*(a0 + r)*(-2*MBH + r)*
     (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)) + 
   (I*\[Omega]*KK[r])/(1 - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r) - 
   ((16*I)*cs^2*Pi*(a0 + r)^2*\[Delta]\[Rho][r])/
    (E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], 
          a0 - M + r]))*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*\[Omega]) + 
   Derivative[1][H1][r], 
  ((-2 + (1 - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r)^(-1))*H0[r])/r + 
   (((-I)*(1 + n))/(r^2*\[Omega]) - 
     (I*E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*r*\[Omega])/
      (2*MBH - r))*H1[r] + ((a0^2*(-3*MBH + r) + 2*a0*r*(-3*MBH + r) + 
      r^2*(-3*MBH + r) - 3*M*(-2*MBH + r)^2)*KK[r])/
    (r*(-2*MBH + r)*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)) - 
   ((4*I)*r*(a0 + r)^4*W[r])/(E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
       (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*
     (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*(-3*M*(-2*MBH + r)^2 + 
      a0^2*(-3*MBH + 2*r) + r^2*(-3*MBH + 2*r) + a0*(-6*MBH*r + 4*r^2))*
     \[Omega]) + Derivative[1][H0][r], 
  -(H0[r]/r) - (I*(1 + n)*H1[r])/(r^2*\[Omega]) + 
   ((a0^2*(-3*MBH + r) + 2*a0*r*(-3*MBH + r) + r^2*(-3*MBH + r) - 
      3*M*(-2*MBH + r)^2)*KK[r])/(r*(-2*MBH + r)*(a0^2 + 4*M*MBH + 2*a0*r - 
      2*M*r + r^2)) - ((4*I)*r*(a0 + r)^4*W[r])/
    (E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], 
          a0 - M + r]))*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*
     (-3*M*(-2*MBH + r)^2 + a0^2*(-3*MBH + 2*r) + r^2*(-3*MBH + 2*r) + 
      a0*(-6*MBH*r + 4*r^2))*\[Omega]) + Derivative[1][KK][r], 
  ((-1/2*I)*M*(a0 + 2*MBH)*(a0^2*(3*MBH - 2*r) + 2*a0*(3*MBH - 2*r)*r + 
      (3*MBH - 2*r)*r^2 + 3*M*(-2*MBH + r)^2)*
     (2*E^(2*Sqrt[M/(2*a0 - M + 4*MBH)]*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], 
          a0 - M + r])*(1 + n)*(2*MBH - r) - 
      E^(Sqrt[M/(2*a0 - M + 4*MBH)]*Pi)*r^3*\[Omega]^2)*H0[r])/
    (E^(2*Sqrt[M/(2*a0 - M + 4*MBH)]*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], 
        a0 - M + r])*r^5*(a0 + r)^5*\[Omega]) + 
   ((I/2)*E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
       (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*M*
     (a0 + 2*MBH)*(a0^2*(3*MBH - 2*r) + 2*a0*(3*MBH - 2*r)*r + 
       (3*MBH - 2*r)*r^2 + 3*M*(-2*MBH + r)^2)^2*\[Omega]*KK[r])/
    ((2*MBH - r)*r^2*(a0 + r)^5*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)) - 
   ((MBH - 2*r + (M*(-2*MBH + r)^2)/(a0 + r)^2 - 
      ((a0^2*(9*MBH - 5*r) + 2*a0*(9*MBH - 5*r)*r + (9*MBH - 5*r)*r^2 + 
         9*M*(-2*MBH + r)^2)*(a0^3*MBH + 3*a0^2*MBH*r + MBH*r^3 + 
         M*r*(12*MBH^2 - 8*MBH*r + r^2) + a0*(4*M*MBH^2 - M*r^2 + 
           3*MBH*r^2)))/((a0 + r)^3*(-3*M*(-2*MBH + r)^2 + 
         a0^2*(-3*MBH + 2*r) + r^2*(-3*MBH + 2*r) + a0*(-6*MBH*r + 4*r^2))))*
     W[r])/(r^2*(1 - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r)) - 
   ((2*I)*Pi*(a0^2*(3*MBH - 2*r) + 2*a0*(3*MBH - 2*r)*r + (3*MBH - 2*r)*r^2 + 
      3*M*(-2*MBH + r)^2)*(2*cs^2*E^(2*Sqrt[M/(2*a0 - M + 4*MBH)]*
         ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])*(1 + n)*
       (2*MBH - r) + E^(Sqrt[M/(2*a0 - M + 4*MBH)]*Pi)*r^3*\[Omega]^2)*
     \[Delta]\[Rho][r])/(E^(2*Sqrt[M/(2*a0 - M + 4*MBH)]*
       ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])*(2*MBH - r)*r^3*
     (a0 + r)^2*\[Omega]) + Derivative[1][W][r], 
  (M*(a0 + 2*MBH)*(a0^2*(3*MBH - r) + 2*a0*(3*MBH - r)*r + (3*MBH - r)*r^2 + 
      3*M*(-2*MBH + r)^2)*H0[r])/(4*csr^2*Pi*r^3*(a0 + r)^3*
     (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)) + 
   ((I/4)*E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
       (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*M*
     (a0 + 2*MBH)*(((1 + n)*(a0^2*(MBH - r) + 2*a0*(MBH - r)*r + 
         (MBH - r)*r^2 + M*(-2*MBH + r)^2))/E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])) - 
      r^3*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*\[Omega]^2)*H1[r])/
    (csr^2*Pi*r^4*(a0 + r)^3*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*
     \[Omega]) - (M*(a0 + 2*MBH)*(a0^2*(MBH - r) + 2*a0*(MBH - r)*r + 
      (MBH - r)*r^2 + M*(-2*MBH + r)^2)*(a0^2*(3*MBH - r) + 
      2*a0*(3*MBH - r)*r + (3*MBH - r)*r^2 + 3*M*(-2*MBH + r)^2)*KK[r])/
    (4*csr^2*Pi*(2*MBH - r)*r^3*(a0 + r)^3*
     (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)^2) - 
   ((I/4)*(-(((2*MBH - r)*r^2*(a0 + r)^5)/E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
          (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))) + 
      ((2*MBH - r)^3*(a0 + r)*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)^2)/
       E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])) - 
      (4*(2*MBH - r)*(a0^2*(MBH - r) + 2*a0*(MBH - r)*r + (MBH - r)*r^2 + 
         M*(-2*MBH + r)^2)*(a0^3*MBH + 3*a0^2*MBH*r + MBH*r^3 + 
         M*r*(12*MBH^2 - 8*MBH*r + r^2) + a0*(4*M*MBH^2 - M*r^2 + 
           3*MBH*r^2)))/E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])) + 
      4*r^4*(a0 + r)^3*(-2*MBH + r)*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*
       \[Omega]^2)*W[r])/(csr^2*Pi*r^3*(a0 + r)^5*
     (r - 2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))^2*
     (4 - (6*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r)*\[Omega]) + 
   ((1 + csr^2 + (-1 - 4*cs^2 + 3*csr^2)*
       (1 - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r))*\[Delta]\[Rho][r])/
    (2*csr^2*(r - 2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))) + 
   Derivative[1][\[Delta]\[Rho]][r]}, 
 {((-I)*ang*Sqrt[Pi/2]*
    Sqrt[1/(E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*
       (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2))]*
    (8*r^5*(a0 + r)^5*\[ScriptM]*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
       Pi/2, 0] + (ang*(-(((2*MBH - r)*r^2*(a0 + r)^5)/
          E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 2*ArcTan[Sqrt[
                M*(2*a0 - M + 4*MBH)], a0 - M + r]))) + 
        (4*(2*MBH - r)*(a0^2*MBH + 2*a0*MBH*r + MBH*r^2 + M*(-2*MBH + r)^2)*
          (a0^3*MBH + 3*a0^2*MBH*r + MBH*r^3 + M*r*(12*MBH^2 - 8*MBH*r + 
             r^2) + a0*(4*M*MBH^2 - M*r^2 + 3*MBH*r^2)))/
         E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 
            2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])) + 
        (a0 + r)*(-2*MBH + r)*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*
         (((2*MBH - r)*(2*M*(-2*MBH + r)^2 + a0^2*(2*MBH + r) + 
             2*a0*r*(2*MBH + r) + r^2*(2*MBH + r)))/
           E^(Sqrt[M/(2*a0 - M + 4*MBH)]*(Pi - 2*ArcTan[Sqrt[
                 M*(2*a0 - M + 4*MBH)], a0 - M + r])) - 4*r^4*(a0 + r)^2*
           \[Omega]^2))*((-1 + \[ScriptM])*\[ScriptM]*SphericalHarmonicY[
          \[ScriptL], \[ScriptM], Pi/2, 0] + 
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(1 + \[ScriptL] + \[ScriptM])*(2 + \[ScriptL] + \[ScriptM])*
           Gamma[1 + \[ScriptL] - \[ScriptM]]]*SphericalHarmonicY[\[ScriptL], 
          2 + \[ScriptM], Pi/2, 0]))/
      (enr*n*(1 - (2*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r)*\[Omega])))/
   (E^((r - rpart)^2/(2*\[Sigma]^2))*(1 + n)*r^7*(a0 + r)^4*\[Sigma]), 
  (2*ang^2*Sqrt[2*Pi]*
    Sqrt[1/(E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*
       (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2))]*(a0^2*MBH + 2*a0*MBH*r + 
     MBH*r^2 + M*(-2*MBH + r)^2)*((-1 + \[ScriptM])*\[ScriptM]*
      SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
     Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
      Sqrt[(2 + 3*\[ScriptL] + \[ScriptL]^2 + 3*\[ScriptM] + 
         2*\[ScriptL]*\[ScriptM] + \[ScriptM]^2)*
        Gamma[1 + \[ScriptL] - \[ScriptM]]]*SphericalHarmonicY[\[ScriptL], 
       2 + \[ScriptM], Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*
    (1 + n)*r^4*(a0 + r)*\[Sigma]), 
  (-2*ang^2*Sqrt[2*Pi]*(-2*MBH + r)*
    Sqrt[1/(E^(Sqrt[-(M/(-2*a0 + M - 4*MBH))]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*
       (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2))]*(a0^2 + 4*M*MBH + 2*a0*r - 
     2*M*r + r^2)*((-1 + \[ScriptM])*\[ScriptM]*SphericalHarmonicY[
       \[ScriptL], \[ScriptM], Pi/2, 0] + 
     Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
      Sqrt[(2 + 3*\[ScriptL] + \[ScriptL]^2 + 3*\[ScriptM] + 
         2*\[ScriptL]*\[ScriptM] + \[ScriptM]^2)*
        Gamma[1 + \[ScriptL] - \[ScriptM]]]*SphericalHarmonicY[\[ScriptL], 
       2 + \[ScriptM], Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
    (n + n^2)*r^4*(a0 + r)*\[Sigma]), 
  (I*ang^2*E^(-1/2*(r - rpart)^2/\[Sigma]^2 + Sqrt[M/(2*a0 - M + 4*MBH)]*
       (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*M*
    (a0 + 2*MBH)*Sqrt[Pi/2]*(2*MBH - r)*
    Sqrt[1/(E^(Sqrt[-(M/(-2*a0 + M - 4*MBH))]*
         (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*
       (a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2))]*
    (4 - (6*(MBH + (M*(-2*MBH + r)^2)/(a0 + r)^2))/r)*
    (((1 + n)*(2*MBH - r)*(a0^2*MBH + 2*a0*MBH*r + MBH*r^2 + 
        M*(-2*MBH + r)^2))/E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
        (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])) + 
     r^3*(-2*MBH + r)*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*\[Omega]^2)*
    ((-1 + \[ScriptM])*\[ScriptM]*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
       Pi/2, 0] + Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
      Sqrt[(4 + 2*n + 2*\[ScriptL]*(1 + \[ScriptM]) + 
         \[ScriptM]*(3 + \[ScriptM]))*Gamma[1 + \[ScriptL] - \[ScriptM]]]*
      SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
   (enr*(n + n^2)*(1 - (2*MBH)/r)*r^8*(a0 + r)^4*\[Sigma]*\[Omega]), 0}, 
 (-2*\[Omega]*(-I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])*Kh[0])/((-2 + r)*
    (1 + n - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])) + 
  (\[Omega]*
    (4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*Sqrt[M/(4 + 2*a0 - M)]*
      (-4 - 2*a0 + M)*\[Omega]*
      (1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]) + Sqrt[(4 + 2*a0 - M)*M]*
      (4*a0*(1 + 2*n - (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) + a0^2*(1 + 2*n - 
         (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
       4*(I + E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(12 + M)*\[Omega] + 
         (4*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*(-8 + M)*\[Omega]^2 + 
         n*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))))*Kh[0])/
   ((2 + a0)^2*Sqrt[(4 + 2*a0 - M)*M]*(I + I*n + 
     2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
    (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])) + (-2 + r)*
   (((2*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M - (10*I)*(2 + a0)*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
        (4 + 2*a0 - M)*M*\[Omega] - I*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
        (4 + 2*a0 - M)*(2*a0*(1 + n)*(12 + 5*M + 18*n) + 
         6*a0^2*(2 + 5*n + 3*n^2) + a0^3*(2 + 5*n + 3*n^2) + 
         4*(4 + 10*n + 6*n^2 + M*(6 + 5*n)))*\[Omega] + 
       24*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M*\[Omega]^2 - 
       E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*
        (6*a0^2*(15 + 36*n + 16*n^2) + a0^3*(15 + 36*n + 16*n^2) + 
         4*a0*(45 + 108*n + 48*n^2 + 2*M*(13 + 9*n)) + 
         8*(15 + 36*n + 16*n^2 + M*(29 + 18*n)))*\[Omega]^2 - 
       (160*I)*(2 + a0)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M*\[Omega]^3 + (8*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*(56 + 76*M + 80*n + 40*M*n + 16*n^2 + 
         6*a0^2*(7 + 10*n + 2*n^2) + a0^3*(7 + 10*n + 2*n^2) + 
         4*a0*(21 + 30*n + 6*n^2 + 5*M*(2 + n)))*\[Omega]^3 - 
       128*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
            6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M*\[Omega]^4 + 
       16*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*
        (8*a0*M*(2 + n) + 12*a0*(9 + 4*n) + 6*a0^2*(9 + 4*n) + 
         a0^3*(9 + 4*n) + 8*(9 + M + 4*n + 2*M*n))*\[Omega]^4 - 
       (128*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
           4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (2*a0^4 + a0*(64 - 20*M) - a0^3*(-16 + M) - 6*a0^2*(-8 + M) + 
         4*(8 - 6*M + M^2))*\[Omega]^5)*Kh[0])/
     (2*(2 + a0)^3*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4 + 2*a0 - M)*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*
      (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2) + 
    (ang*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 - (Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(-8*enr*n*(I + E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]) + ang*(-1 + \[ScriptM])*\[Omega]*
          (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*(3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (8*enr*n*(1 + n)*\[Sigma]*
      (-1 + (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) + (-2 + r)^2*
   (((-6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(72*Sqrt[(4 + 2*a0 - M)*M^3] - 
         (2 + a0)*Sqrt[(4 + 2*a0 - M)*M]*(178 + 52*n + a0*(35 + 26*n))) + 
       (6*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4 + 2*a0 - M)*M^2*
        (540*Sqrt[(4 + 2*a0 - M)*M^3] - Sqrt[(4 + 2*a0 - M)*M]*
          (a0^2*(545 + 374*n) + 4*a0*(761 + 374*n) + 4*(977 - 9*M + 374*n)))*
        \[Omega] + (3*I)*(2 + a0)*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
        ((4 + 2*a0 - M)*M)^(3/2)*(3*(2 + a0)^5*n*(19 + 29*n + 10*n^2) + 
         32*M^2*(71 + 62*n + 31*a0*(1 + n)) + 2*(2 + a0)*M*
          (-56 + 796*n + 744*n^2 + a0^2*(151 + 337*n + 186*n^2) + 
           2*a0*(191 + 536*n + 372*n^2)))*\[Omega] + 
       24*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(67*Sqrt[(4 + 2*a0 - M)*M^3] - 
         Sqrt[(4 + 2*a0 - M)*M]*(3280 - 21*M + 1984*n + a0^2*(754 + 496*n) + 
           4*a0*(787 + 496*n)))*\[Omega]^2 + 3*(2 + a0)*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*((4 + 2*a0 - M)*M)^(3/2)*
        (10*a0^4*(-171 + 202*n + 476*n^2 + 184*n^3) + 
         a0^5*(-171 + 202*n + 476*n^2 + 184*n^3) + 
         32*(-171 + 202*n + 476*n^2 + 184*n^3 + M^2*(874 + 604*n) + 
           M*(-423 + 654*n + 736*n^2)) + 16*a0*(M^2*(796 + 604*n) + 
           24*M*(4 + 117*n + 92*n^2) + 5*(-171 + 202*n + 476*n^2 + 
             184*n^3)) + 8*a0^3*(M*(354 + 750*n + 368*n^2) + 
           5*(-171 + 202*n + 476*n^2 + 184*n^3)) + 
         8*a0^2*(3*M*(409 + 1218*n + 736*n^2) + 10*(-171 + 202*n + 476*n^2 + 
             184*n^3)))*\[Omega]^2 + (192*I)*(2 + a0)^2*
        E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*M^2*
        (225*Sqrt[(4 + 2*a0 - M)*M^3] + Sqrt[(4 + 2*a0 - M)*M]*
          (172 + 15*M + 584*n + a0^2*(223 + 146*n) + a0*(532 + 584*n)))*
        \[Omega]^3 - (12*I)*(2 + a0)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        ((4 + 2*a0 - M)*M)^(3/2)*(30*a0^4*(-73 - 26*n + 44*n^2 + 24*n^3) + 
         3*a0^5*(-73 - 26*n + 44*n^2 + 24*n^3) + 
         32*(M^2*(934 + 544*n) + M*(-899 + 214*n + 480*n^2) + 
           3*(-73 - 26*n + 44*n^2 + 24*n^3)) + 
         24*a0^2*(3*M*(35 + 286*n + 160*n^2) + 10*(-73 - 26*n + 44*n^2 + 
             24*n^3)) + 16*a0*(16*M^2*(57 + 34*n) + 
           6*M*(-193 + 268*n + 240*n^2) + 15*(-73 - 26*n + 44*n^2 + 
             24*n^3)) + 8*a0^3*(M*(287 + 590*n + 240*n^2) + 
           15*(-73 - 26*n + 44*n^2 + 24*n^3)))*\[Omega]^3 + 
       1536*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(87*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(-164 + 5*M + 72*n + 2*a0^2*(14 + 9*n) + 
           a0*(-26 + 72*n)))*\[Omega]^4 + 
       16*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*M*
        (54*Sqrt[(4 + 2*a0 - M)*M^7] - Sqrt[4 + 2*a0 - M]*M^(5/2)*
          (8708 - 196*M + 5400*n + 675*a0^2*(3 + 2*n) + a0*(8404 + 5400*n)) - 
         3*Sqrt[4 + 2*a0 - M]*M^(3/2)*(478*M^2 + 8*M*(2283 + 1310*n) + 
           48*(69 + 220*n + 100*n^2) + 24*a0^3*(213 + 308*n + 100*n^2) + 
           a0^4*(783 + 1012*n + 300*n^2) + 16*a0*(702 + 1496*n + 600*n^2 + 
             M*(1184 + 655*n)) + 2*a0^2*(M*(2453 + 1310*n) + 
             12*(495 + 836*n + 300*n^2))) - Sqrt[-(M*(-4 - 2*a0 + M))]*
          (36*a0^5*(-45 - 58*n + 4*n^2 + 8*n^3) + 
           3*a0^6*(-45 - 58*n + 4*n^2 + 8*n^3) - 
           4*(296*M^3 - M^2*(3133 + 1926*n) + M*(34980 + 12432*n - 2544*
                n^2) - 48*(-45 - 58*n + 4*n^2 + 8*n^3)) + 
           48*a0^3*(M*(-418 + 115*n + 106*n^2) + 10*(-45 - 58*n + 4*n^2 + 8*
                n^3)) + 3*a0^4*(M*(97 + 652*n + 212*n^2) + 
             60*(-45 - 58*n + 4*n^2 + 8*n^3)) + 4*a0*(M^2*(4139 + 1926*n) + 
             48*M*(-1171 - 307*n + 106*n^2) + 144*(-45 - 58*n + 4*n^2 + 8*
                n^3)) + 3*a0^2*(M^2*(1715 + 642*n) + 24*M*(-1649 - 192*n + 
               212*n^2) + 240*(-45 - 58*n + 4*n^2 + 8*n^3))))*\[Omega]^4 - 
       (1536*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(90*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(-220 + 6*M + 24*n + 4*a0*(-19 + 6*n) + 
           a0^2*(17 + 6*n)))*\[Omega]^5 + 
       (192*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
           4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*
        M*(15*Sqrt[(4 + 2*a0 - M)*M^7] + Sqrt[4 + 2*a0 - M]*M^(5/2)*
          (2188 + 175*M + 968*n + 4*a0*(557 + 242*n) + a0^2*(567 + 242*n)) + 
         Sqrt[4 + 2*a0 - M]*M^(3/2)*(405*M^2 + 24*M*(135 + 94*n) + 
           3*a0^4*(93 + 76*n + 12*n^2) + 16*a0^3*(84 + 95*n + 18*n^2) + 
           16*(-105 + 76*n + 36*n^2) + 16*a0*(-48 + 295*M + 228*n + 141*M*n + 
             72*n^2) + 2*a0^2*(804 + 775*M + 1824*n + 282*M*n + 432*n^2)) - 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(595*M^3 + 192*(-1 + 2*n) + 
           36*a0^5*(-1 + 2*n) + a0^6*(-3 + 6*n) + 4*M^2*(653 + 38*n) - 
           16*M*(-651 - 316*n + 12*n^2) + 4*a0*(144*(-1 + 2*n) + 
             M^2*(585 + 38*n) - 24*M*(-217 - 76*n + 4*n^2)) - 
           3*a0^4*(60 - 120*n + M*(-13 + 12*n + 4*n^2)) + 
           a0^2*(720*(-1 + 2*n) + M^2*(517 + 38*n) - 24*M*(-549 - 140*n + 12*
                n^2)) - 8*a0^3*(60 - 120*n + M*(-345 - 52*n + 12*n^2))))*
        \[Omega]^5 - 6144*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(7*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(-16 - 4*a0 + 2*a0^2 + M))*\[Omega]^6 + 
       12288*E^((Sqrt[M/(4 + 2*a0 - M)]*(7*Pi - 
            10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(a0^4*Sqrt[-(M*(-4 - 2*a0 + M))]*(2 + n) + 
         a0^3*Sqrt[-(M*(-4 - 2*a0 + M))]*(-23 + 2*n) - 
         8*(8*Sqrt[(4 + 2*a0 - M)*M] + 5*Sqrt[(4 + 2*a0 - M)*M^3] + 
           4*Sqrt[(4 + 2*a0 - M)*M]*n - 2*Sqrt[(4 + 2*a0 - M)*M^3]*n) + 
         4*a0^2*Sqrt[-(M*(-4 - 2*a0 + M))]*(M*(2 + n) - 3*(11 + n)) + 
         4*a0*Sqrt[-(M*(-4 - 2*a0 + M))]*(-47 - 10*n + M*(-1 + 4*n)))*
        \[Omega]^6 + (49152*I)*(2 + a0)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (2*Pi - 3*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)^(3/2)*M^(5/2)*(6*a0 + 3*a0^2 + 4*M)*\[Omega]^7)*Kh[0])/
     (36*(2 + a0)^6*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*((4 + 2*a0 - M)*M)^(3/2)*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(-3 + (16*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega] + 16*E^(Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         \[Omega]^2)^2) + 
    (ang*Pi*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*
          (3 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
          ((10*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
             Sqrt[4 + 2*a0 - M]*M^(3/2)*Sqrt[2/Pi]*\[Omega]^2*
             (-2 + (7*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                     (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
              4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/\[Sigma] + 
           Sqrt[(4 + 2*a0 - M)*M]*((2*M*Sqrt[2/Pi]*(3 - (8*I)*
                 E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                10*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                (13*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                       (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                 \[Omega]^3 + 12*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                 \[Omega]^4))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
             (2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2*
              ((-4*Sqrt[2/Pi]*(2 - rpart)*(-5 + (14*I)*
                   E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                  8*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/
                (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
               (9 - (76*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] - 32*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + n*(-10 + (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]))/(E^((2 - rpart)^2/
                   (2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma])))) + 
         8*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*enr*n*\[Omega]*
          ((2*M*Sqrt[M/(4 + 2*a0 - M)]*(-4 - 2*a0 + M)*Sqrt[2/Pi]*
             (-9 + (45*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                     (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
              59*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 - 
              (20*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3))/
            (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
           Sqrt[(4 + 2*a0 - M)*M]*((-2*(2 + a0)^2*Sqrt[2/Pi]*(2 - rpart)*(
                -9 + (36*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                       (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                 \[Omega] + 44*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                 \[Omega]^2 - (16*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]^3))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^
                3) + (68 - (468*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega] - 560*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + (128*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]^3 + 8*n*(-5 + (9*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2) + 4*M*(-9 + (27*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 29*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 - (12*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3) + 4*a0*(17 - (117*I)*
                  E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                 140*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                 (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^3 + 2*n*(-5 + (9*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                      (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))*\[Omega]^2)) + a0^2*(17 - (117*I)*
                  E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                 140*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                 (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^3 + 2*n*(-5 + (9*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                      (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))*\[Omega]^2)))/(E^((2 - rpart)^2/
                 (2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]))))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*(3 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        ((10*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
           Sqrt[4 + 2*a0 - M]*M^(3/2)*Sqrt[2/Pi]*\[Omega]^2*
           (-2 + (7*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
            4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/\[Sigma] + 
         Sqrt[(4 + 2*a0 - M)*M]*((2*M*Sqrt[2/Pi]*
             (3 - (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                     (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
              10*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
              (13*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
              12*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^4))/
            (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
           (2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2*
            ((-4*Sqrt[2/Pi]*(2 - rpart)*(-5 + (14*I)*
                 E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                8*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/
              (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
             (9 - (76*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                      (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega] - 32*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + n*(-10 + (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]))/(E^((2 - rpart)^2/
                 (2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]))))*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (48*(2 + a0)^2*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*enr*
      Sqrt[(4 + 2*a0 - M)*M]*n*(1 + n)*\[Omega]*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2)), 
 (((2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega] - 8*E^(Sqrt[M/(4 + 2*a0 - M)]*
        (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2)*
    Kh[0])/((-2 + r)*(1 + n - 
     (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])) + 
  ((4*E^(Sqrt[M/(4 + 2*a0 - M)]*
        (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
      Sqrt[4 + 2*a0 - M]*M^(3/2)*\[Omega]^2*
      (-I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
     E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*Sqrt[(4 + 2*a0 - M)*M]*\[Omega]*
      (4*(1 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (1 + 2*n - (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) + 4*a0*(1 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*(1 + 2*n - (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) + a0^2*(1 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*(1 + 2*n - (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) + M*(4 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega] + 48*E^(Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
          \[Omega]^2)))*Kh[0])/((2 + a0)^2*Sqrt[(4 + 2*a0 - M)*M]*
    (1 + n - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])*
    (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])) + (-2 + r)*
   ((((2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4 + 2*a0 - M)*M*(1 + n)*
        (2*M + (2 + a0)^2*n) - (2*I)*(2 + a0)^2*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*M^2*\[Omega] - 
       I*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M*(6*a0^2*(-1 + 9*n + 7*n^2) + 
         a0^3*(-1 + 9*n + 7*n^2) + 4*(-2 + 9*M + 18*n + 5*M*n + 14*n^2) + 
         2*a0*(-6 + 11*M + 54*n + 5*M*n + 42*n^2))*\[Omega] - 
       8*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*\[Omega]^2 - (2 + a0)*
        E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[(4 + 2*a0 - M)*M]*
        (2*Sqrt[4 + 2*a0 - M]*M^(3/2)*(20 + 8*n + a0*(11 + 4*n)) + 
         Sqrt[(4 + 2*a0 - M)*M]*(6*a0^2*(3 + 36*n + 16*n^2) + 
           a0^3*(3 + 36*n + 16*n^2) + 8*(3 + M + 36*n - 4*M*n + 16*n^2) + 
           a0*(36 + 26*M + 432*n - 16*M*n + 192*n^2)))*\[Omega]^2 - 
       (32*I)*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*\[Omega]^3 + 
       (8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
            6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        ((4*a0^3 - 2*a0^2*(-10 + M) - 8*a0*(-4 + M) + (-4 + M)^2)*M + 
         3*(4 + 2*a0 - M)*M^2 + (4 + 2*a0 - M)*(8*a0^3*(7 + 10*n + 2*n^2) + 
           a0^4*(7 + 10*n + 2*n^2) - 2*(M^2 + 6*M*(7 + 4*n) - 
             8*(7 + 10*n + 2*n^2)) + a0^2*(-2*M*(13 + 6*n) + 
             24*(7 + 10*n + 2*n^2)) + a0*(-2*M*(47 + 24*n) + 
             32*(7 + 10*n + 2*n^2))))*\[Omega]^3 - 128*(2 + a0)^2*
        E^(Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 4*ArcTan[(2 + a0 - M)/
              Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*M^2*\[Omega]^4 + 
       16*(2 + a0)*E^(Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
           4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*
        M*(72 + 32*n - 16*M*(4 + n) - 8*a0*M*(7 + n) + 12*a0*(9 + 4*n) + 
         6*a0^2*(9 + 4*n) + a0^3*(9 + 4*n))*\[Omega]^4 - 
       (128*I)*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        (2*a0^4 + a0^2*(48 - 22*M) - a0^3*(-16 + M) + 8*(4 - 5*M + M^2) + 
         a0*(64 - 60*M + 8*M^2))*\[Omega]^5)*Kh[0])/
     (2*(2 + a0)^4*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4 + 2*a0 - M)*M*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*
      (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2) + 
    (ang*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 - (Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(-8*E^(Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*enr*n*
          \[Omega] + ang*(-1 + \[ScriptM])*
          (1 - I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
           4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*(1 - I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
         4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2)*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (8*enr*n*(1 + n)*\[Sigma]*
      (-1 + (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) + (-2 + r)^2*
   ((((-36*I)*(2 + a0)^4*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        (-(M*(-4 - 2*a0 + M)))^(3/2)*(1 + n) - 240*(2 + a0)^4*
        E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*M*(-(M*(-4 - 2*a0 + M)))^(3/2)*(1 + n)*
        \[Omega] - 9*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
        ((4 + 2*a0 - M)*M)^(3/2)*(1 + n)*(40*M^2 + 3*(2 + a0)^4*n*
          (-1 + 2*n) + 2*(2 + a0)*M*(-2 + 24*n + a0*(5 + 12*n)))*\[Omega] - 
       (6*I)*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(24*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(a0^2*(5 + 38*n) + 4*a0*(23 + 38*n) + 
           4*(41 + 9*M + 38*n)))*\[Omega]^2 + (3*I)*(2 + a0)*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*M*
        (108*(2 + a0)*Sqrt[4 + 2*a0 - M]*M^(3/2)*(1 + n)*
          (2*M + (2 + a0)^2*n) + Sqrt[(4 + 2*a0 - M)*M]*
          (90*a0^4*(3 - 21*n - n^2 + 14*n^3) + 9*a0^5*(3 - 21*n - n^2 + 
             14*n^3) + 16*(2*M^2*(62 + 53*n) + M*(-143 + 241*n + 234*n^2) + 
             18*(3 - 21*n - n^2 + 14*n^3)) + 36*a0^2*
            (M*(-3 + 115*n + 78*n^2) + 20*(3 - 21*n - n^2 + 14*n^3)) + 
           8*a0*(2*M^2*(71 + 53*n) + 3*M*(-85 + 293*n + 234*n^2) + 
             90*(3 - 21*n - n^2 + 14*n^3)) + 2*a0^3*
            (M*(85 + 397*n + 234*n^2) + 180*(3 - 21*n - n^2 + 14*n^3))))*
        \[Omega]^2 - 24*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(35*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(976 + 57*M + 832*n + 2*a0^2*(89 + 104*n) + 
           a0*(844 + 832*n)))*\[Omega]^3 + 3*(2 + a0)*
        E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))*((4 + 2*a0 - M)*M)^(3/2)*
        (10*a0^4*(405 - 902*n - 244*n^2 + 280*n^3) + 
         a0^5*(405 - 902*n - 244*n^2 + 280*n^3) + 
         32*(405 - 902*n - 244*n^2 + 280*n^3 + 2*M^2*(104 + 77*n) + 
           M*(-519 + 444*n + 376*n^2)) + 16*a0*(22*M^2*(13 + 7*n) + 
           3*M*(-463 + 480*n + 376*n^2) + 5*(405 - 902*n - 244*n^2 + 
             280*n^3)) + 8*a0^2*(3*M*(-329 + 516*n + 376*n^2) + 
           10*(405 - 902*n - 244*n^2 + 280*n^3)) + 
         4*a0^3*(M*(-117 + 552*n + 376*n^2) + 10*(405 - 902*n - 244*n^2 + 
             280*n^3)))*\[Omega]^3 + (192*I)*(2 + a0)^2*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*M^2*
        (-9*Sqrt[(4 + 2*a0 - M)*M^3] + Sqrt[(4 + 2*a0 - M)*M]*
          (412 - 3*M + 296*n + 4*a0*(103 + 74*n) + a0^2*(103 + 74*n)))*
        \[Omega]^4 - (12*I)*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        ((4 + 2*a0 - M)*M)^(3/2)*(30*a0^4*(131 - 114*n - 20*n^2 + 24*n^3) + 
         a0^5*(393 - 342*n - 60*n^2 + 72*n^3) - 32*(-393 + 342*n + 60*n^2 - 
           72*n^3 + 4*M^2*(1 + 8*n) + M*(221 - 118*n - 96*n^2)) - 
         16*a0*(2*M^2*(-9 + 16*n) - 6*M*(-141 + 28*n + 48*n^2) + 
           15*(-131 + 114*n + 20*n^2 - 24*n^3)) + 
         24*a0^2*(M*(-321 - 6*n + 96*n^2) + 10*(131 - 114*n - 20*n^2 + 
             24*n^3)) + 8*a0^3*(M*(-169 - 34*n + 48*n^2) + 
           15*(131 - 114*n - 20*n^2 + 24*n^3)))*\[Omega]^4 + 
       1536*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(-9*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(64 - 13*M + 48*n + a0^2*(31 + 12*n) + 
           a0*(94 + 48*n)))*\[Omega]^5 - 
       16*E^(Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
           4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*
        M*(270*Sqrt[(4 + 2*a0 - M)*M^7] - Sqrt[4 + 2*a0 - M]*M^(5/2)*
          (9*a0^2*(27 + 98*n) + 4*a0*(491 + 882*n) + 
           4*(739 - 86*M + 882*n)) - 3*Sqrt[4 + 2*a0 - M]*M^(3/2)*
          (a0^4*(747 + 164*n - 84*n^2) - 24*a0^3*(-235 - 48*n + 28*n^2) + 
           16*a0*(1242 + 208*n - 168*n^2 + M*(130 + 53*n)) + 
           2*(4632 + 91*M^2 + 672*n - 672*n^2 + M*(996 + 424*n)) + 
           2*a0^2*(M*(271 + 106*n) + 12*(663 + 124*n - 84*n^2))) + 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(36*a0^5*(99 - 58*n + 4*n^2 + 8*n^3) + 
           3*a0^6*(99 - 58*n + 4*n^2 + 8*n^3) - 3*a0^2*(5*M^2*(125 + 78*n) + 
             24*M*(-433 + 264*n + 52*n^2) - 240*(99 - 58*n + 4*n^2 + 8*
                n^3)) - 4*a0*(M^2*(253 + 1170*n) + 48*M*(-356 + 55*n + 26*
                n^2) - 144*(99 - 58*n + 4*n^2 + 8*n^3)) - 
           3*a0^4*(M*(485 + 572*n + 52*n^2) - 60*(99 - 58*n + 4*n^2 + 8*
                n^3)) - 48*a0^3*(M*(-17 + 209*n + 26*n^2) - 
             10*(99 - 58*n + 4*n^2 + 8*n^3)) + 
           4*(-17*M^3 + M^2*(1369 - 1170*n) - 12*M*(-871 - 44*n + 52*n^2) + 
             48*(99 - 58*n + 4*n^2 + 8*n^3))))*\[Omega]^5 - 
       (1536*I)*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(-6*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(-28 - 30*M + 24*n + 4*a0*(11 + 6*n) + 
           a0^2*(29 + 6*n)))*\[Omega]^6 + 
       (192*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(7*Pi - 
            10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M*(33*Sqrt[(4 + 2*a0 - M)*M^7] - Sqrt[4 + 2*a0 - M]*
          M^(5/2)*(3348 - 57*M + 1144*n + 11*a0^2*(75 + 26*n) + 
           4*a0*(831 + 286*n)) + Sqrt[4 + 2*a0 - M]*M^(3/2)*
          (-816 + 328*M + 123*M^2 + 64*n - 48*M*n + 192*n^2 - 
           48*a0*(60 + M*(-4 + n) + 4*n - 8*n^2) + 
           3*a0^4*(-89 - 12*n + 4*n^2) + 16*a0^3*(-99 - 13*n + 6*n^2) - 
           2*a0^2*(1668 + 192*n - 144*n^2 + M*(-7 + 6*n))) - 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(213*M^3 + 192*(-1 + 2*n) + 
           36*a0^5*(-1 + 2*n) + a0^6*(-3 + 6*n) - 4*M^2*(1139 + 10*n) + 
           16*M*(-207 + 52*n + 12*n^2) - 4*a0*(144 - 288*n + 
             M^2*(719 + 10*n) - 24*M*(-63 + 32*n + 4*n^2)) + 
           3*a0^4*(60*(-1 + 2*n) + M*(159 + 76*n + 4*n^2)) + 
           a0^2*(720*(-1 + 2*n) - M^2*(299 + 10*n) + 24*M*(-69 + 140*n + 12*
                n^2)) + 8*a0^3*(60*(-1 + 2*n) + M*(153 + 184*n + 12*n^2))))*
        \[Omega]^6 + 6144*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (2*Pi - 3*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(-Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(16 + 4*a0 - 2*a0^2 + 5*M))*\[Omega]^7 - 
       256*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(2*Pi - 
           3*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*
        M^2*(12*a0^4*Sqrt[-(M*(-4 - 2*a0 + M))]*(55 + 8*n) + 
         24*a0^3*Sqrt[-(M*(-4 - 2*a0 + M))]*(109 + 26*n) + 
         384*(-Sqrt[(4 + 2*a0 - M)*M] - 7*Sqrt[(4 + 2*a0 - M)*M^3] + 
           Sqrt[(4 + 2*a0 - M)*M]*n + Sqrt[(4 + 2*a0 - M)*M^3]*n) + 
         48*a0^2*Sqrt[-(M*(-4 - 2*a0 + M))]*(51 + 30*n + 2*M*(2 + n)) + 
         96*a0*Sqrt[-(M*(-4 - 2*a0 + M))]*(-5 + 14*n + 2*M*(-5 + 2*n)))*
        \[Omega]^7 + (24576*I)*(2 + a0)*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 14*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(6*a0^2 + 3*a0^3 - 8*M)*
        (4 + 2*a0 - M)^(3/2)*M^(5/2)*\[Omega]^8)*Kh[0])/
     (36*(2 + a0)^6*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*((4 + 2*a0 - M)*M)^(3/2)*
      \[Omega]*(I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(-3 + (16*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega] + 16*E^(Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         \[Omega]^2)^2) + 
    (ang*Pi*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
          (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
          ((2*M*Sqrt[M/(4 + 2*a0 - M)]*(-4 - 2*a0 + M)*Sqrt[2/Pi]*
             (-3 + (11*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                     (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
              19*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 - 
              (4*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3))/
            (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
           Sqrt[(4 + 2*a0 - M)*M]*(((2*I)*(2 + a0)^2*Sqrt[2/Pi]*(2 - rpart)*
               (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2*(-3*I + 
                4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))/
              (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
             (4*a0*(-27 + (75*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 68*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3 + 2*n*(-3 + (5*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                   4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2)) + a0^2*
                (-27 + (75*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 68*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3 + 2*n*(-3 + (5*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                   4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2)) + 4*
                (-27 + (75*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 68*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3 + 2*n*(-3 + (5*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                   4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
                 M*(-2 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                    \[Omega] - 5*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2 + (12*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3)))/(E^((2 - rpart)^2/
                 (2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]))) - 
         8*enr*n*((2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[
                 M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                     -(M*(-4 - 2*a0 + M))]]))*M*Sqrt[M/(4 + 2*a0 - M)]*
             (-4 - 2*a0 + M)*Sqrt[2/Pi]*\[Omega]^2*(12*I + 
              19*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
              (4*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/\[Sigma] + 
           Sqrt[(4 + 2*a0 - M)*M]*((-4*(2 + a0)^2*E^(-1/2*(2 - rpart)^2/
                   \[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[
                2/Pi]*(2 - rpart)*\[Omega]^2*(3*I + 
                10*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/
              \[Sigma]^3 + (4*a0*(-6*I - 16*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (39*I)*
                  E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                 4*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
                 (32*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                       (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^4 + 2*n*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (9*I)*E^(Sqrt[M/(4 + 2*a0 - 
                         M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 
                          2*a0 + M))]]))*\[Omega]^2 + 
                   4*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3)) + 
               a0^2*(-6*I - 16*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + (39*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + 4*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^3 + (32*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^4 + 2*n*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (9*I)*E^(Sqrt[M/(4 + 2*a0 - 
                         M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 
                          2*a0 + M))]]))*\[Omega]^2 + 
                   4*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3)) + 4*
                (-6*I - 16*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] - (3*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  (-13 + M)*\[Omega]^2 - E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*(-4 + 13*M)*\[Omega]^3 + (4*I)*
                  E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))*(8 + 3*M)*\[Omega]^4 + 
                 2*n*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                    \[Omega] + (9*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2 + 4*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                        2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                      2)*\[Omega]^3)))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[
                2*Pi]*\[Sigma]))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
         Pi/2, 0] + ang*\[Omega]*(3*I + 
         4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        ((2*M*Sqrt[M/(4 + 2*a0 - M)]*(-4 - 2*a0 + M)*Sqrt[2/Pi]*
           (-3 + (11*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                   (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
            19*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 - 
            (4*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3))/
          (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
         Sqrt[(4 + 2*a0 - M)*M]*(((2*I)*(2 + a0)^2*Sqrt[2/Pi]*(2 - rpart)*
             (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2*
             (-3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))/
            (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
           (4*a0*(-27 + (75*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega] + 68*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]^3 + 2*n*(-3 + (5*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2)) + a0^2*(-27 + (75*I)*
                E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 68*
                E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + (32*I)*
                E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 2*n*
                (-3 + (5*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2)) + 4*(-27 + (75*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega] + 68*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + (32*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]^3 + 2*n*(-3 + (5*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2) + M*(-2 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] - 5*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (12*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3)))/
            (E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma])))*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (48*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*enr*
      Sqrt[(4 + 2*a0 - M)*M]*n*(1 + n)*\[Omega]*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2)), 
 Kh[0] + 
  ((-2 + r)*
    (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*Sqrt[4 + 2*a0 - M]*M^(3/2)*
      \[Omega]*(I + 18*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega] - (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        \[Omega]^2 + n*(I + 12*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])) - Sqrt[(4 + 2*a0 - M)*M]*
      ((2 + a0)^2*(n + n^2 + (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*(1 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) + 2*M*(1 + n - (7*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega] - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*n*
          \[Omega] + 10*E^(Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
          \[Omega]^2 + 12*E^(Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*n*
          \[Omega]^2 - (40*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]^3)))*Kh[0])/((2 + a0)^2*Sqrt[(4 + 2*a0 - M)*M]*
    (1 + n - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])*
    (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
           2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^
     2) + (-2 + r)^2*
   (((-2*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        (1 + n) + (10*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*M*(1 + n)*
        \[Omega] + I*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*(1 + n)*
        (32*M^2 + 3*(2 + a0)^4*(-1 + n)*n + 2*(2 + a0)*M*(-4 + a0 + 18*n + 
           9*a0*n))*\[Omega] - 4*(2 + a0)^2*
        E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*M*(5 + 6*n)*\[Omega]^2 + 
       E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*
        (8*a0^3*(9 - 43*n - 18*n^2 + 16*n^3) + a0^4*(9 - 43*n - 18*n^2 + 
           16*n^3) + 16*(9 - 43*n - 18*n^2 + 16*n^3 + 2*M^2*(11 + 8*n) + 
           M*(-26 + 24*n + 30*n^2)) + 8*a0*(M*(-31 + 63*n + 60*n^2) + 
           4*(9 - 43*n - 18*n^2 + 16*n^3)) + 4*a0^2*(M*(-5 + 39*n + 30*n^2) + 
           6*(9 - 43*n - 18*n^2 + 16*n^3)))*\[Omega]^2 + 
       (16*I)*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*M*
        (9 + 10*n)*\[Omega]^3 - (2*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (8*a0^3*(51 - 80*n - 48*n^2 + 8*n^3) + a0^4*(51 - 80*n - 48*n^2 + 
           8*n^3) + 8*a0^2*(-1 + 2*n)*(7*M*(2 + n) + 
           3*(-51 - 22*n + 4*n^2)) + 16*(51 - 80*n - 48*n^2 + 8*n^3 + 
           M^2*(38 + 20*n) + 2*M*(-38 + 9*n + 14*n^2)) + 
         32*a0*(51 - 80*n - 48*n^2 + 8*n^3 + M*(-26 + 15*n + 14*n^2)))*
        \[Omega]^3 + 64*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        (3 + 2*n)*\[Omega]^4 + 
       16*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(96*a0*(-7 + 3*n + 2*n^2) + 
         24*a0^3*(-7 + 3*n + 2*n^2) + 3*a0^4*(-7 + 3*n + 2*n^2) - 
         4*a0*M*(-43 - 2*n + 8*n^2) - 4*a0^2*(M*(-7 + n + 2*n^2) - 
           18*(-7 + 3*n + 2*n^2)) - 8*(2*M^2*(5 + 2*n) - 
           6*(-7 + 3*n + 2*n^2) + M*(-29 - 4*n + 4*n^2)))*\[Omega]^4 - 
       (256*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*M*
        \[Omega]^5 + (32*I)*(2 + a0)*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (72 + 54*a0^2 + 9*a0^3 - 16*M*(-3 + n) + a0*(108 - 8*M*n))*
        \[Omega]^5 + 3072*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*M*
        \[Omega]^6)*Kh[0])/(4*(2 + a0)^4*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
      \[Omega]*(I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*
      (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2*
      (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])) + (ang*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 - 
        (Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
              Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(3*I + I*n + 
           2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
          (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
         8*enr*n*(1 + n - I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega] - I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*n*
            \[Omega] - 2*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                 (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*(3*I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*(3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega])*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (16*enr*n*(1 + n)*\[Sigma]*\[Omega]*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) + (-2 + r)^3*
   (((6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(1 + n)*(36*Sqrt[(4 + 2*a0 - M)*M^3] - 
         Sqrt[(4 + 2*a0 - M)*M]*(316 - 18*M + 52*n + 13*a0^2*(4 + n) + 
           a0*(262 + 52*n))) - (6*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*
        (4 + 2*a0 - M)*M^2*(1 + n)*(306*Sqrt[(4 + 2*a0 - M)*M^3] - 
         Sqrt[(4 + 2*a0 - M)*M]*(3036 - 96*M + 748*n + a0^2*(543 + 187*n) + 
           4*a0*(651 + 187*n)))*\[Omega] - 
       (3*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*((4 + 2*a0 - M)*M)^(3/2)*(1 + n)*
        (2304*M^3 + 15*(2 + a0)^6*(-2 + n)^2*n + 4*(2 + a0)*M^2*
          (-556 + 626*n + a0*(55 + 313*n)) + 2*(2 + a0)^2*M*
          (8*(-5 - 128*n + 60*n^2) + a0^2*(-4 - 79*n + 120*n^2) + 
           2*a0*(-68 - 335*n + 240*n^2)))*\[Omega] - 
       6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(2*Sqrt[4 + 2*a0 - M]*M^(3/2)*(169 + 181*n) - 
         Sqrt[(4 + 2*a0 - M)*M]*(282*M*(1 + n) + a0^2*(1693 + 2664*n + 
             992*n^2) + 4*a0*(1777 + 2730*n + 992*n^2) + 
           4*(1861 + 2796*n + 992*n^2)))*\[Omega]^2 + 
       3*E^((Sqrt[M/(4 + 2*a0 - M)]*(3*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*((4 + 2*a0 - M)*M)^(3/2)*
        (a0^6*(180 - 893*n - 60*n^2 + 516*n^3 - 92*n^4) - 
         12*a0^5*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4) - 
         32*(24*M^3*(40 + 31*n) + M^2*(-2299 + 15*n + 1396*n^2) + 
           2*M*(257 - 1744*n - 889*n^2 + 458*n^3) + 
           2*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4)) - 
         16*a0*(M^2*(-2273 + 1851*n + 2792*n^2) + 
           2*M*(281 - 6463*n - 2719*n^2 + 1832*n^3) + 
           12*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4)) - 
         4*a0^4*(M*(-40 - 889*n - 52*n^2 + 458*n^3) + 
           15*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4)) - 
         16*a0^2*(M^2*(13 + 918*n + 698*n^2) + 3*M*(-83 - 2861*n - 941*n^2 + 
             916*n^3) + 15*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4)) - 
         8*a0^3*(M*(-313 - 4753*n - 1045*n^2 + 1832*n^3) + 
           20*(-180 + 893*n + 60*n^2 - 516*n^3 + 92*n^4)))*\[Omega]^2 - 
       (24*I)*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (Pi - ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(Sqrt[4 + 2*a0 - M]*M^(3/2)*(965 + 948*n) + 
         Sqrt[(4 + 2*a0 - M)*M]*(M*(471 + 444*n) + 
           4*a0^2*(29 + 160*n + 146*n^2) + 8*(-269 - 40*n + 292*n^2) + 
           4*a0*(-211 + 280*n + 584*n^2)))*\[Omega]^3 - 
       (3*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - ArcTan[(2 + a0 - M)/
             Sqrt[-(M*(-4 - 2*a0 + M))]]))*((4 + 2*a0 - M)*M)^(3/2)*
        (a0^6*(2139 - 5236*n - 1088*n^2 + 1904*n^3 - 144*n^4) - 
         12*a0^5*(-2139 + 5236*n + 1088*n^2 - 1904*n^3 + 144*n^4) - 
         64*(-2139 + 5236*n + 1088*n^2 - 1904*n^3 + 144*n^4 + 
           180*M^3*(13 + 8*n) + M^2*(-6774 - 550*n + 2216*n^2) + 
           2*M*(1546 - 4021*n - 2070*n^2 + 564*n^3)) - 
         64*a0*(-6417 + 15708*n + 3264*n^2 - 5712*n^3 + 432*n^4 + 
           M^2*(-3757 + 1330*n + 2216*n^2) + M*(3211 - 15544*n - 6852*n^2 + 
             2256*n^3)) - 16*a0^3*(M*(49 - 12832*n - 3996*n^2 + 2256*n^3) + 
           10*(-2139 + 5236*n + 1088*n^2 - 1904*n^3 + 144*n^4)) - 
         4*a0^4*(2*M*(-35 - 2665*n - 642*n^2 + 564*n^3) + 
           15*(-2139 + 5236*n + 1088*n^2 - 1904*n^3 + 144*n^4)) - 
         16*a0^2*(M^2*(-740 + 3210*n + 2216*n^2) + 
           6*M*(583 - 7230*n - 2712*n^2 + 1128*n^3) + 
           15*(-2139 + 5236*n + 1088*n^2 - 1904*n^3 + 144*n^4)))*\[Omega]^3 - 
       192*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
            6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(3*Sqrt[4 + 2*a0 - M]*M^(3/2)*(145 + 152*n) + 
         Sqrt[(4 + 2*a0 - M)*M]*(M*(29 + 8*n) + 4*(-445 - 396*n + 72*n^2) + 
           4*a0*(-307 - 258*n + 72*n^2) + a0^2*(-169 - 120*n + 72*n^2)))*
        \[Omega]^4 + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (5*Pi - 6*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M*(540*Sqrt[4 + 2*a0 - M]*M^(7/2)*(1 + n) - 
         24*Sqrt[4 + 2*a0 - M]*M^(3/2)*(3*a0^3 - 6*a0^2*(-2 + M) + 
           4*(11 - 6*M)*M + 2*a0*(6 + 5*M))*(193 + 136*n) - 
         4*Sqrt[4 + 2*a0 - M]*M^(5/2)*(5162 - 7756*n - 8748*n^2 - 
           2*M*(8297 + 5972*n) + a0*(259 - 11384*n - 8748*n^2) - 
           27*a0^2*(43 + 139*n + 81*n^2)) + 3*Sqrt[4 + 2*a0 - M]*M^(3/2)*
          (4*M^2*(8981 + 3569*n) + 8*M*(-14985 + 2232*n + 6268*n^2) + 
           8*a0*M*(-11217 + 4136*n + 6268*n^2) + 2*a0^2*M*(-7449 + 6040*n + 
             6268*n^2) + 96*a0*(1083 - 1958*n - 1208*n^2 + 312*n^3) + 
           3*a0^4*(403 - 1586*n - 684*n^2 + 312*n^3) + 
           48*a0^2*(1279 - 2751*n - 1550*n^2 + 468*n^3) + 
           16*a0^3*(939 - 2565*n - 1288*n^2 + 468*n^3) + 
           16*(3951 - 6246*n - 4148*n^2 + 936*n^3)) + 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(36*a0^5*(-2559 + 3532*n + 1248*n^2 - 
             624*n^3 + 16*n^4) + 3*a0^6*(-2559 + 3532*n + 1248*n^2 - 
             624*n^3 + 16*n^4) - 4*a0*(M^2*(49427 + 17480*n + 7956*n^2) - 
             48*M*(2280 - 4549*n - 1868*n^2 + 92*n^3) - 
             144*(-2559 + 3532*n + 1248*n^2 - 624*n^3 + 16*n^4)) - 
           6*a0^2*(M^2*(5061 + 2686*n + 1326*n^2) - 24*M*(987 - 4465*n - 1578*
                n^2 + 92*n^3) - 120*(-2559 + 3532*n + 1248*n^2 - 624*n^3 + 16*
                n^4)) - 8*(M^3*(4772 + 5345*n) + M^2*(34244 + 9422*n + 3978*
                n^2) - 6*M*(7989 - 8930*n - 4316*n^2 + 184*n^3) - 
             24*(-2559 + 3532*n + 1248*n^2 - 624*n^3 + 16*n^4)) + 
           24*a0^3*(M*(231 - 8426*n - 2576*n^2 + 184*n^3) + 
             20*(-2559 + 3532*n + 1248*n^2 - 624*n^3 + 16*n^4)) + 
           3*a0^4*(M*(-669 - 7586*n - 1996*n^2 + 184*n^3) + 
             60*(-2559 + 3532*n + 1248*n^2 - 624*n^3 + 16*n^4))))*
        \[Omega]^4 + (1536*I)*(2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(3*Sqrt[4 + 2*a0 - M]*M^(3/2)*(21 + 23*n) + 
         Sqrt[(4 + 2*a0 - M)*M]*(-(M*(19 + 15*n)) + 4*(-82 - 62*n + 3*n^2) + 
           a0^2*(-31 - 26*n + 3*n^2) + 2*a0*(-113 - 88*n + 6*n^2)))*
        \[Omega]^5 + (16*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (3*Pi - 4*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M*(96*Sqrt[4 + 2*a0 - M]*M^(3/2)*
          (3*a0^3 - 6*a0^2*(-2 + M) + 4*(11 - 6*M)*M + 2*a0*(6 + 5*M))*
          (36 + 19*n) - 18*Sqrt[4 + 2*a0 - M]*M^(7/2)*(296 + 191*n) + 
         Sqrt[4 + 2*a0 - M]*M^(5/2)*(59740 + 5184*n - 19632*n^2 - 
           2*M*(29425 + 17643*n) - 3*a0^2*(-777 + 2176*n + 1636*n^2) - 
           4*a0*(-8633 + 2616*n + 4908*n^2)) + Sqrt[-(M*(-4 - 2*a0 + M))]*
          (36*a0^5*(1497 - 1004*n - 448*n^2 + 64*n^3) + 
           3*a0^6*(1497 - 1004*n - 448*n^2 + 64*n^3) + 
           48*a0^3*(M*(523 + 309*n - 6*n^2 + 12*n^3) + 
             10*(1497 - 1004*n - 448*n^2 + 64*n^3)) + 
           3*a0^4*(M*(1475 + 770*n - 20*n^2 + 24*n^3) + 
             60*(1497 - 1004*n - 448*n^2 + 64*n^3)) + 
           2*(M^3*(25723 + 12855*n) + M^2*(-63050 - 20256*n + 8520*n^2) + 
             24*M*(263 + 234*n + 12*n^2 + 24*n^3) + 96*(1497 - 1004*n - 448*
                n^2 + 64*n^3)) + 4*a0*(M^2*(-16667 - 1320*n + 4260*n^2) + 
             48*M*(220 + 175*n + 2*n^2 + 12*n^3) + 144*(1497 - 1004*n - 448*
                n^2 + 64*n^3)) + 3*a0^2*(M^2*(-603 + 2496*n + 1420*n^2) + 
             24*M*(701 + 478*n - 4*n^2 + 24*n^3) + 240*(1497 - 1004*n - 448*
                n^2 + 64*n^3))) - 3*Sqrt[4 + 2*a0 - M]*M^(3/2)*
          (16*a0^3*(1140 - 995*n - 686*n^2 + 60*n^3) + 
           a0^4*(1647 - 1958*n - 1188*n^2 + 120*n^3) - 
           2*(M^2*(5578 + 4857*n) - 4*M*(-6861 - 592*n + 948*n^2) + 
             8*(-4299 + 1966*n + 1924*n^2 - 120*n^3)) + 
           8*a0*(M*(-6619 - 648*n + 948*n^2) + 8*(1803 - 997*n - 870*n^2 + 60*
                n^3)) + 2*a0^2*(M*(-6377 - 704*n + 948*n^2) + 
             12*(2933 - 2002*n - 1556*n^2 + 120*n^3))))*\[Omega]^5 + 
       1536*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)*M^2*(2*Sqrt[4 + 2*a0 - M]*M^(3/2)*(22 + 13*n) - 
         Sqrt[(4 + 2*a0 - M)*M]*(300 + 80*n + 2*M*(26 + 5*n) + 
           a0^2*(27 + 8*n) + 4*a0*(51 + 14*n)))*\[Omega]^6 + 
       64*E^((Sqrt[M/(4 + 2*a0 - M)]*(7*Pi - 10*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*(4 + 2*a0 - M)*M*
        (48*Sqrt[4 + 2*a0 - M]*M^(3/2)*(3*a0^3 - 6*a0^2*(-2 + M) + 
           4*(11 - 6*M)*M + 2*a0*(6 + 5*M))*(29 + 10*n) - 
         3*Sqrt[4 + 2*a0 - M]*M^(7/2)*(1493 + 902*n) + Sqrt[4 + 2*a0 - M]*
          M^(5/2)*(40724 + 15680*n - 2448*n^2 - M*(15739 + 6946*n) + 
           a0*(28564 + 9136*n - 2448*n^2) + a0^2*(4101 + 648*n - 612*n^2)) + 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(M^3*(9991 + 3322*n) - 
           576*(-131 + 32*n + 16*n^2) - 108*a0^5*(-131 + 32*n + 16*n^2) - 
           9*a0^6*(-131 + 32*n + 16*n^2) - 48*M*(-2171 + 194*n + 92*n^2) + 
           4*M^2*(-19271 - 5888*n + 204*n^2) - 
           3*a0^4*(180*(-131 + 32*n + 16*n^2) + M*(-911 + 74*n + 52*n^2)) - 
           24*a0^3*(60*(-131 + 32*n + 16*n^2) + M*(-857 + 158*n + 62*n^2)) - 
           3*a0^2*(M^2*(1165 + 184*n - 68*n^2) + 720*(-131 + 32*n + 16*n^2) + 
             24*M*(-1049 + 206*n + 72*n^2)) + 
           4*a0*(-432*(-131 + 32*n + 16*n^2) - 24*M*(-1487 + 218*n + 82*
                n^2) + M^2*(-11383 - 3220*n + 204*n^2))) + 
         3*Sqrt[4 + 2*a0 - M]*M^(3/2)*(M^2*(12625 + 5182*n) + 
           8*M*(-125 - 152*n + 4*n^2) + a0^4*(-723 + 194*n + 164*n^2) + 
           8*a0^3*(-1005 + 146*n + 186*n^2) + 16*(-1863 - 70*n + 252*n^2) + 
           8*a0*(-6300 + 56*n + 920*n^2 + M*(839 + 204*n + 4*n^2)) + 
           2*a0^2*(M*(1803 + 560*n + 4*n^2) + 12*(-1289 + 86*n + 208*n^2))))*
        \[Omega]^6 + (6144*I)*(2 + a0)^2*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
          (2*Pi - 3*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        (4 + 2*a0 - M)*M^2*(-7*Sqrt[(4 + 2*a0 - M)*M^3] + 
         Sqrt[(4 + 2*a0 - M)*M]*(40 + 28*a0 + 4*a0^2 + 11*M))*\[Omega]^7 - 
       (256*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(2*Pi - 
           3*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*(4 + 2*a0 - M)*
        M*(48*Sqrt[4 + 2*a0 - M]*M^(3/2)*(3*a0^3 - 6*a0^2*(-2 + M) + 
           4*(11 - 6*M)*M + 2*a0*(6 + 5*M))*(6 + n) - 3*Sqrt[4 + 2*a0 - M]*
          M^(7/2)*(103 + 102*n) + Sqrt[4 + 2*a0 - M]*M^(5/2)*
          (7804 - 1181*M + 3264*n + 6*M*n + 144*n^2 + 
           3*a0^2*(349 + 132*n + 12*n^2) + 4*a0*(1499 + 606*n + 36*n^2)) - 
         3*Sqrt[4 + 2*a0 - M]*M^(3/2)*(a0^4*(93 + 8*n) + 
           12*a0^3*(97 + 12*n) + 48*(93 + 20*n) - 17*M^2*(163 + 30*n) + 
           8*M*(541 + 92*n) + 8*a0*(954 + 176*n + M*(151 + 18*n)) + 
           a0^2*(48*(97 + 15*n) - 2*M*(239 + 56*n))) + 
         Sqrt[-(M*(-4 - 2*a0 + M))]*(6912 + 1296*a0^5 + 108*a0^6 + 
           a0^3*(17280 + M*(3588 - 384*n)) + a0^4*(6480 + M*(435 - 24*n)) + 
           M^3*(89 - 78*n) - 48*M*(-667 + 4*n) - 4*M^2*(2353 + 696*n + 
             36*n^2) - 3*a0^2*(-8640 + 48*M*(-119 + 9*n) + 
             M^2*(187 + 116*n + 12*n^2)) - 4*a0*(-5184 + 12*M*(-821 + 28*n) + 
             M^2*(1457 + 522*n + 36*n^2))))*\[Omega]^7 - 
       98304*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 
            14*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        (4 + 2*a0 - M)^(3/2)*(12 + 6*a0 - M)*M^(5/2)*\[Omega]^8)*Kh[0])/
     (36*(2 + a0)^6*E^(Sqrt[M/(4 + 2*a0 - M)]*Pi)*((4 + 2*a0 - M)*M)^(3/2)*
      \[Omega]*(I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(I + I*n + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(5*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])*(-3 + (16*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega] + 16*E^(Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         \[Omega]^2)^2) + 
    (ang*Pi*\[Mu]*(\[ScriptM]*(-(ang*(-1 + \[ScriptM])*
           (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
           ((E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*M*
              Sqrt[M/(4 + 2*a0 - M)]*(-4 - 2*a0 + M)*Sqrt[2/Pi]*\[Omega]^2*
              (-40*I - 187*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                      (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega] + (186*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + 40*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega]^3 + I*n*(-10 + (53*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 44*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2)))/\[Sigma] + Sqrt[(4 + 2*a0 - M)*M]*
             ((M*Sqrt[2/Pi]*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + (5*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 - 99*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^3 + (178*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^4 + 72*E^((5*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^5 + n*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (19*I)*E^(Sqrt[M/(4 + 2*a0 - 
                         M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 
                          2*a0 + M))]]))*\[Omega]^2 + 
                   7*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
                   (12*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                         (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^4)))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
                \[Sigma]) + (2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
               \[Omega]^2*((-2*Sqrt[2/Pi]*(2 - rpart)*(4 + n - (2*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
                  (5*I + 14*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                    \[Omega] - (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
                  \[Sigma]^3) + (-180*I - 629*E^((Sqrt[M/(4 + 2*a0 - M)]*
                      (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (548*I)*
                   E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                  128*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - 
                          M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
                  n^2*(5*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                     \[Omega]) + (2*I)*n*(-17 + (66*I)*E^((Sqrt[M/(4 + 2*a0 - 
                          M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 
                          2*a0 + M))]]))/2)*\[Omega] + 40*
                     E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/
                 (E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]))))) + 
         8*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*enr*n*\[Omega]*
          ((Sqrt[4 + 2*a0 - M]*M^(3/2)*Sqrt[2/Pi]*(3*I + 
              4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
             (3*I + 17*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
              (3*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 - 
              10*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                      Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
              n*(3*I + 17*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                       (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                 \[Omega] - (11*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                 \[Omega]^2)))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
           Sqrt[(4 + 2*a0 - M)*M]*(-(((2 + a0)^2*Sqrt[2/Pi]*(2 - rpart)*
                (-3 + (10*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 8*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2)*(3 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] - 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + n*(3 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega])))/(E^((2 - rpart)^2/
                  (2*\[Sigma]^2))*\[Sigma]^3)) + 
             (4*a0*(43 - (191*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                      2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] - 265*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (228*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3 + 128*E^(2*Sqrt[M/(4 + 2*a0 - 
                       M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))*\[Omega]^4 + n^2*(-5 + (9*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                   4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
                 2*n*(19 - (91*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                        2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                      2)*\[Omega] - 114*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2 + (40*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3)) + a0^2*(43 - (191*I)*
                  E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                         Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                 265*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                 (228*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega]^3 + 128*E^(2*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^4 + n^2*(-5 + (9*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + 4*E^(Sqrt[M/(4 + 2*a0 - M)]*
                      (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))*\[Omega]^2) + 2*n*(19 - (91*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                   114*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
                   (40*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                    \[Omega]^3)) + 2*(86 - (382*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] - 530*E^(Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2 + (456*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3 + 256*E^(2*Sqrt[M/(4 + 2*a0 - 
                       M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))*\[Omega]^4 + 2*n^2*(-5 + (9*I)*
                    E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
                   4*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
                 4*n*(19 - (91*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                        2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                      2)*\[Omega] - 114*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2 + (40*I)*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
                       (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]^3) + M*(3*I + 
                   4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                          Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
                  (6*I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                          (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                    \[Omega] - I*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                       2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                    \[Omega]^2 + 18*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                        2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                      2)*\[Omega]^3 + n*(6*I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*
                         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega] + (3*I)*E^(Sqrt[M/(4 + 2*a0 - 
                          M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 
                          2*a0 + M))]]))*\[Omega]^2))))/(E^((2 - rpart)^2/
                 (2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]))))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
       ang*(3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        ((E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*M*
           Sqrt[M/(4 + 2*a0 - M)]*(-4 - 2*a0 + M)*Sqrt[2/Pi]*\[Omega]^2*
           (-40*I - 187*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                   (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
            (186*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2 + 
            40*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]^3 + 
            I*n*(-10 + (53*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega] + 44*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
               \[Omega]^2)))/\[Sigma] + Sqrt[(4 + 2*a0 - M)*M]*
          ((M*Sqrt[2/Pi]*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega] + (5*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
               \[Omega]^2 - 99*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega]^3 + (178*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
               \[Omega]^4 + 72*E^((5*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega]^5 + n*(-3*I - 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega] + (19*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                 \[Omega]^2 + 7*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                 \[Omega]^3 + (12*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                 \[Omega]^4)))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) + 
           (2 + a0)^2*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2*
            ((-2*Sqrt[2/Pi]*(2 - rpart)*(4 + n - (2*I)*
                 E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*(5*I + 
                14*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                        Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] - 
                (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                       Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2))/
              (E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
             (-180*I - 629*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                      (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega] + (548*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                \[Omega]^2 + 128*E^((3*Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                \[Omega]^3 + n^2*(5*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                     (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + 
                          M))]]))/2)*\[Omega]) + (2*I)*n*(-17 + 
                 (66*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[
                        (2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                  \[Omega] + 40*E^(Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                     2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
                  \[Omega]^2))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
               \[Sigma]))))*Conjugate[
         1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (48*(2 + a0)^2*E^((3*Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*enr*
      Sqrt[(4 + 2*a0 - M)*M]*n*(1 + n)*\[Omega]^2*
      (I + 2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[
                -(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
      (3*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])^2*
      (5*I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))), 0, 0, (-I)*K\[Infinity]0*r*\[Omega] + 
  K\[Infinity]0*(n - (2*I)*(1 + M)*\[Omega]) + 
  ((K\[Infinity]0*((3*I)*n^2 + 3*n*(I + 4*(1 + M)*\[Omega] + 
         I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(9 - (24*I)*\[Omega] + 
         2*M^2*\[Omega]*(-12*I + 5*(2 + a0)*\[Omega]) + 
         3*M*(3 + (2*I)*(-2 + 3*a0)*\[Omega] - 2*(-2 + a0 + a0^2)*
            \[Omega]^2))))/(6*\[Omega]) + 
    (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
          \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) + 
    (I/2)*\[Omega]*
     (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)) + 
      (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
           21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
              \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
              \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
             (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
             2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
           15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
             2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
             2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
             a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^
                 3)))))/(30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
        ((-I)*\[ScriptM]*(2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - 
             (9*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*
              \[Omega]^3) + ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
             54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
           \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
          (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
           (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
          Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
          Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
             \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
          SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
       (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
        (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
      (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
             18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
            \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*
                M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                   2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
           0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
              \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
               \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
             \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
            (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
           2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
        n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2))))/r + 
  (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
          21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
           \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
          27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
             \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
          3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
             \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
        I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
          2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*\[Omega] + 
            (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
          5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - 
            (15*I)*(4 + 20*a0 + 9*a0^2)*\[Omega]^2 + 
            2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*\[Omega]^3) + 
          15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 8*a0^4*\[Omega]^3 + 
            2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*
             (-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
            a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
      (\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*
        (2*enr*n*(-36*I - (9*I)*n + 18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*
            \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
         ang*(-1 + \[ScriptM])*\[Omega]*((9*I)*n + \[Omega]*
            (-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*
                \[Omega] + 2*a0*(2 + a0)*\[Omega]^2))))*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
       ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
           3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^2)))*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
      (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
    (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 
            5*n^3*(-I + 8*(1 + M)*\[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
            10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^
                3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + 
                (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + (142 - 85*a0 + 
                  66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - 
              I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
               \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                 \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                 \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
              5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
                 \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*
                 \[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*
                 \[Omega]^4))))/(120*\[Omega]^3) + 
         (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(
                36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                  16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
              2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
                  (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
                  2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
                  M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 
              81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^
                2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
              n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
             Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
             Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
             SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
          (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
           \[Omega]^2)) + \[Omega]*
        ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
            10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
              I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
            5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
              (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
              M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                 \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
            I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
              2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
              5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                  19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                 \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
          (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
           (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                   \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*((3*I)*n^2 + 
                n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
                \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                  M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + 
              n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
              \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                   \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                  \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
                \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - 
                 \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], 
              Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*
           \[Sigma]*\[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
      (2 + a0)*M*\[Omega]^2))/r^2 + 
  ((K\[Infinity]0*((-15*I)*n^3*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
       15*n^2*(-6 + (36*I)*\[Omega] - I*M^2*\[Omega]*
          (-36 - (16*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^2) + 
         M*(-6 - I*(-10 + 31*a0)*\[Omega] + 6*(-4 + a0^2)*\[Omega]^2)) + 
       5*n*\[Omega]*(18*M*(3 + (16*I)*\[Omega])*\[Omega] + 
         54*(3*I + 8*\[Omega]) + 6*a0^3*M*\[Omega]^2*(9*I + 2*M*\[Omega]) + 
         8*M^3*\[Omega]*(54 + (69*I)*\[Omega] - 10*\[Omega]^2) - 
         2*a0^2*M*\[Omega]*(-117 + (18*I)*(1 + 8*M)*\[Omega] + 
           2*M*(-9 + 5*M)*\[Omega]^2) - 6*M^2*(-27*I + 43*\[Omega] - 
           (80*I)*\[Omega]^2 + 8*\[Omega]^3) - 
         a0*M*(162*I + 3*(51 + 259*M)*\[Omega] - (12*I)*(-12 - 28*M + 23*M^2)*
            \[Omega]^2 + 80*M^2*\[Omega]^3)) + 
       \[Omega]^2*(540*(3 - (8*I)*\[Omega]) + 4*M^4*\[Omega]*
          (-1080*I + 852*(2 + a0)*\[Omega] + (65*I)*(2 + a0)^2*\[Omega]^2) + 
         M^3*(1620 + (360*I)*(14 + 31*a0)*\[Omega] + 
           (4188 - 9972*a0 - 6033*a0^2)*\[Omega]^2 - (5*I)*(2 + a0)^2*
            (-22 + 69*a0)*\[Omega]^3) - 90*M*(4*a0^4*\[Omega]^2 - 
           4*a0^3*\[Omega]*(3*I + \[Omega]) + a0^2*(-12 + (15*I)*\[Omega] - 
             8*\[Omega]^2) - 4*(2 - (3*I)*\[Omega] + 8*\[Omega]^2) + 
           a0*(-1 - (12*I)*\[Omega] + 16*\[Omega]^2)) + 
         15*M^2*(16*a0^3*(12 + I*\[Omega])*\[Omega]^2 + 
           (6*I)*a0^4*\[Omega]^3 + a0^2*\[Omega]*(-502*I + 129*\[Omega]) - 
           2*a0*(93 + (56*I)*\[Omega] + 190*\[Omega]^2) + 
           4*(-12 + (14*I)*\[Omega] + 65*\[Omega]^2 + (8*I)*\[Omega]^3)))))/
     (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (I*ang*Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(2*enr*n*((-6*I)*(4 + n + (2*I)*\[Omega]) + 
           M^2*\[Omega]*(12 - (8*I)*(2 + a0)*\[Omega] + (2 + a0)^2*
              \[Omega]^2) + M*(-24*I - (6*I)*n - 2*(-10 + a0)*\[Omega] + 
             (2 + a0)*n*\[Omega] + (2*I)*(-4 + 4*a0 + 3*a0^2)*\[Omega]^2)) - 
         ang*(-1 + \[ScriptM])*\[Omega]*
          (n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
           \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - (4*I)*
                \[Omega] + (3*I)*a0^2*\[Omega]) + M^2*(12 - (8*I)*(2 + a0)*
                \[Omega] + (2 + a0)^2*\[Omega]^2))))*SphericalHarmonicY[
         \[ScriptL], \[ScriptM], Pi/2, 0] - ang*\[Omega]*
        (n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
         \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - (4*I)*\[Omega] + 
             (3*I)*a0^2*\[Omega]) + M^2*(12 - (8*I)*(2 + a0)*\[Omega] + 
             (2 + a0)^2*\[Omega]^2)))*Conjugate[
         1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
      (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (\[Omega]*(-((2 + a0)*M*\[Omega]*
         ((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + I*
                (2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
               I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - (2 + a0)*
                (-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 5*(1 + M)*
                \[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - 
               (192 + (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 
                 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 
                 52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
             \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 2*
                M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
                 (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 15*(-12 + 
                 (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 5*
                M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - (4*I)*
                  (54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 
                   123*a0^2 + 78*a0^3)*\[Omega]^3) - 5*M*(36 - 
                 (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
                  \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*
                  \[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*
                  \[Omega]^4))))/(120*\[Omega]^3) + 
          (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                   16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                  \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                      \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                   (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                 \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                     3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                     2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*SphericalHarmonicY[
               \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
              (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^2))) + 
       (3*I)*((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
            10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
              I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
            5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
              (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
              M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                 \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
            I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
              2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
              5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                  19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                 \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
          (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
           (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                   \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*((3*I)*n^2 + 
                n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
                \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                  M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + 
              n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
              \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                   \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                  \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
                \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - 
                 \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], 
              Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*
           \[Sigma]*\[Omega]^3) + 2*(1 + M)*
          ((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
                I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*
                 \[Omega] - I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 
                5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
                 \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
                  9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
                 (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
                 \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                  77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                  (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                   \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                  (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                   \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                  2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
                5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 
                    90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 
                    8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 
                    2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
           (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                 (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 
                    3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 
                    3*M)*M*\[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                       \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                    (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                  \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                      3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                      2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*
               SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
              ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
                I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*
                 (-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
                 (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
                1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*Sqrt[
                (2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + \[ScriptL]*
                   (3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
               SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
            (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
             \[Omega]^2)))))/(9 - (6*I)*(1 + M)*\[Omega] + 
      (2 + a0)*M*\[Omega]^2))/r^3, (-I)*K\[Infinity]0*\[Omega] + 
  (E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*K\[Infinity]0*
    r^2*(a0 + r)*\[Omega]^2)/((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + 
      r^2]) + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*K\[Infinity]0*
    r*(a0 + r)*\[Omega]*(n - (2*I)*(1 + M)*\[Omega]))/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  ((K\[Infinity]0*((3*I)*n^2 + 3*n*(I + 4*(1 + M)*\[Omega] + 
         I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(9 - (24*I)*\[Omega] + 
         2*M^2*\[Omega]*(-12*I + 5*(2 + a0)*\[Omega]) + 
         3*M*(3 + (2*I)*(-2 + 3*a0)*\[Omega] - 2*(-2 + a0 + a0^2)*
            \[Omega]^2))))/(6*\[Omega]) + 
    (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
          \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) + 
    (I/2)*\[Omega]*
     (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)) + 
      (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
           21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
              \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
              \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
             (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
             2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
           15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
             2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
             2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
             a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^
                 3)))))/(30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
        ((-I)*\[ScriptM]*(2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - 
             (9*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*
              \[Omega]^3) + ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
             54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
           \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
          (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
           (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
          Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
          Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
             \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
          SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
       (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
        (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
      (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
             18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
            \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*
                M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                   2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
           0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
              \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
               \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
             \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
            (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
           2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
        n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2))))/r^2 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*((3*I)*n^2 + 3*n*(I + 4*(1 + M)*\[Omega] + 
          I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(9 - (24*I)*\[Omega] + 
          2*M^2*\[Omega]*(-12*I + 5*(2 + a0)*\[Omega]) + 
          3*M*(3 + (2*I)*(-2 + 3*a0)*\[Omega] - 2*(-2 + a0 + a0^2)*
             \[Omega]^2))))/(6*\[Omega]) + 
     (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
           \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
        ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) + 
     (I/2)*\[Omega]*(-1/30*(K\[Infinity]0*(45*n^3 - 
           5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 21*(2 + a0)*M*\[Omega]^2 + 
             (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
           5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 27*(3*I + 8*\[Omega]) + 
             M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*\[Omega] + (20 - 92*a0 - 
                 51*a0^2)*\[Omega]^2) + 3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - 
               (2*I)*(-12 + 8*a0 + 7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 
                 3*a0^2)*\[Omega]^3)) + I*\[Omega]^2*
            (-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
             2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*\[Omega] + 
               (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
             5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                 9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
                \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
               8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*
                a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + a0*
                (-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
         (\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
            21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
               \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 
                7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
               \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - (9*I)*(-36 + 28*a0 + 
                23*a0^2)*\[Omega]^2 + 2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*
               \[Omega]^3) + 15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(
                -3*I + 2*\[Omega]) + 2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 
                4*\[Omega]^2) + 2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
              a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
           (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
               \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
            ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
              54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
            \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
           (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
            (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
           Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
           Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
              \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
           SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
        (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
              18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
             \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
                3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                    2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
            0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
               \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                 \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
              \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
             (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
            2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*
         enr*n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)))))/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (2*(-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
           21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
              \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
              \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - 
             (15*I)*(4 + 20*a0 + 9*a0^2)*\[Omega]^2 + 
             2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*\[Omega]^3) + 
           15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 8*a0^4*\[Omega]^3 + 
             2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*
              (-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
             a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
       (\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
            18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
            2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
           \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
              3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                  2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
        ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
            3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^2)))*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
     (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*
                \[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
             10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*
                \[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*
                (8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + 
                 (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*
                \[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - 
                 a0*(24 + 71*M))*\[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*
                (-240*I + 77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*
                (-360 - (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 
                   653*a0^2)*\[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 
                 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*
                (33*I + (22 - 205*a0)*\[Omega] - (4*I)*(54 - 53*a0 + 
                   104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) - 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - 
                 (142 + 35*a0 + 90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 
                   19*a0^2 + 8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + 
                   a0^3 + 2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
          (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                   16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                  \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                      \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                   (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                 \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                     3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                     2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*SphericalHarmonicY[
               \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
              (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^2)) + \[Omega]*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/r^3 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*
    (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
           21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
              \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
              \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - 
             (15*I)*(4 + 20*a0 + 9*a0^2)*\[Omega]^2 + 
             2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*\[Omega]^3) + 
           15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 8*a0^4*\[Omega]^3 + 
             2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*
              (-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
             a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
       (\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
            18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
            2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
           \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
              3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                  2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
        ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
            3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^2)))*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
     (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*
                \[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
             10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*
                \[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*
                (8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + 
                 (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*
                \[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - 
                 a0*(24 + 71*M))*\[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*
                (-240*I + 77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*
                (-360 - (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 
                   653*a0^2)*\[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 
                 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*
                (33*I + (22 - 205*a0)*\[Omega] - (4*I)*(54 - 53*a0 + 
                   104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) - 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - 
                 (142 + 35*a0 + 90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 
                   19*a0^2 + 8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + 
                   a0^3 + 2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
          (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                   16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                  \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                      \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                   (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                 \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                     3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                     2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*SphericalHarmonicY[
               \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
              (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^2)) + \[Omega]*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/((-2 + r)*r*
    Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (3*((K\[Infinity]0*((-15*I)*n^3*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
        15*n^2*(-6 + (36*I)*\[Omega] - I*M^2*\[Omega]*
           (-36 - (16*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^2) + 
          M*(-6 - I*(-10 + 31*a0)*\[Omega] + 6*(-4 + a0^2)*\[Omega]^2)) + 
        5*n*\[Omega]*(18*M*(3 + (16*I)*\[Omega])*\[Omega] + 
          54*(3*I + 8*\[Omega]) + 6*a0^3*M*\[Omega]^2*(9*I + 2*M*\[Omega]) + 
          8*M^3*\[Omega]*(54 + (69*I)*\[Omega] - 10*\[Omega]^2) - 
          2*a0^2*M*\[Omega]*(-117 + (18*I)*(1 + 8*M)*\[Omega] + 
            2*M*(-9 + 5*M)*\[Omega]^2) - 6*M^2*(-27*I + 43*\[Omega] - 
            (80*I)*\[Omega]^2 + 8*\[Omega]^3) - 
          a0*M*(162*I + 3*(51 + 259*M)*\[Omega] - (12*I)*(-12 - 28*M + 
              23*M^2)*\[Omega]^2 + 80*M^2*\[Omega]^3)) + 
        \[Omega]^2*(540*(3 - (8*I)*\[Omega]) + 4*M^4*\[Omega]*
           (-1080*I + 852*(2 + a0)*\[Omega] + (65*I)*(2 + a0)^2*\[Omega]^2) + 
          M^3*(1620 + (360*I)*(14 + 31*a0)*\[Omega] + 
            (4188 - 9972*a0 - 6033*a0^2)*\[Omega]^2 - (5*I)*(2 + a0)^2*
             (-22 + 69*a0)*\[Omega]^3) - 90*M*(4*a0^4*\[Omega]^2 - 
            4*a0^3*\[Omega]*(3*I + \[Omega]) + a0^2*(-12 + (15*I)*\[Omega] - 
              8*\[Omega]^2) - 4*(2 - (3*I)*\[Omega] + 8*\[Omega]^2) + 
            a0*(-1 - (12*I)*\[Omega] + 16*\[Omega]^2)) + 
          15*M^2*(16*a0^3*(12 + I*\[Omega])*\[Omega]^2 + 
            (6*I)*a0^4*\[Omega]^3 + a0^2*\[Omega]*(-502*I + 129*\[Omega]) - 
            2*a0*(93 + (56*I)*\[Omega] + 190*\[Omega]^2) + 
            4*(-12 + (14*I)*\[Omega] + 65*\[Omega]^2 + (8*I)*\[Omega]^3)))))/
      (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (I*ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*((-6*I)*(4 + n + (2*I)*\[Omega]) + 
            M^2*\[Omega]*(12 - (8*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^
                2) + M*(-24*I - (6*I)*n - 2*(-10 + a0)*\[Omega] + 
              (2 + a0)*n*\[Omega] + (2*I)*(-4 + 4*a0 + 3*a0^2)*\[Omega]^2)) - 
          ang*(-1 + \[ScriptM])*\[Omega]*
           (n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
            \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - 
                (4*I)*\[Omega] + (3*I)*a0^2*\[Omega]) + M^2*(12 - 
                (8*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^2))))*
         SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
        ang*\[Omega]*(n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
          \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - (4*I)*\[Omega] + 
              (3*I)*a0^2*\[Omega]) + M^2*(12 - (8*I)*(2 + a0)*\[Omega] + 
              (2 + a0)^2*\[Omega]^2)))*Conjugate[
          1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (\[Omega]*(-((2 + a0)*M*\[Omega]*
          ((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
                I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*
                 \[Omega] - I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 
                5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
                 \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
                  9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
                 (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
                 \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                  77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                  (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                   \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                  (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                   \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                  2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
                5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 
                    90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 
                    8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 
                    2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
           (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                 (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 
                    3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 
                    3*M)*M*\[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                       \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                    (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                  \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                      3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                      2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*
               SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
              ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
                I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*
                 (-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
                 (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
                1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*Sqrt[
                (2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + \[ScriptL]*
                   (3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
               SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
            (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
             \[Omega]^2))) + (3*I)*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3) + 2*(1 + M)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*
                (-I + 8*(1 + M)*\[Omega] + I*(2 + a0)*M*\[Omega]^2) + 10*n^2*
                (5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 24*M^2)*
                  \[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*
                (12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
                  \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
                   9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
                  (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
                  \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                   77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                   (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                    \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                   (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                    \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                   2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
                 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 
                     90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 
                     8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 
                     2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
            (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                  (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 
                     3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 
                     3*M)*M*\[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                        \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                     (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                   \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                       3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                       2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*
                SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + ang*
                \[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
                 I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*
                  (-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
                  (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
                 1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
                Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                   \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
                SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
             (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
              \[Omega]^2)))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/r^4 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*
    ((K\[Infinity]0*((-15*I)*n^3*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
        15*n^2*(-6 + (36*I)*\[Omega] - I*M^2*\[Omega]*
           (-36 - (16*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^2) + 
          M*(-6 - I*(-10 + 31*a0)*\[Omega] + 6*(-4 + a0^2)*\[Omega]^2)) + 
        5*n*\[Omega]*(18*M*(3 + (16*I)*\[Omega])*\[Omega] + 
          54*(3*I + 8*\[Omega]) + 6*a0^3*M*\[Omega]^2*(9*I + 2*M*\[Omega]) + 
          8*M^3*\[Omega]*(54 + (69*I)*\[Omega] - 10*\[Omega]^2) - 
          2*a0^2*M*\[Omega]*(-117 + (18*I)*(1 + 8*M)*\[Omega] + 
            2*M*(-9 + 5*M)*\[Omega]^2) - 6*M^2*(-27*I + 43*\[Omega] - 
            (80*I)*\[Omega]^2 + 8*\[Omega]^3) - 
          a0*M*(162*I + 3*(51 + 259*M)*\[Omega] - (12*I)*(-12 - 28*M + 
              23*M^2)*\[Omega]^2 + 80*M^2*\[Omega]^3)) + 
        \[Omega]^2*(540*(3 - (8*I)*\[Omega]) + 4*M^4*\[Omega]*
           (-1080*I + 852*(2 + a0)*\[Omega] + (65*I)*(2 + a0)^2*\[Omega]^2) + 
          M^3*(1620 + (360*I)*(14 + 31*a0)*\[Omega] + 
            (4188 - 9972*a0 - 6033*a0^2)*\[Omega]^2 - (5*I)*(2 + a0)^2*
             (-22 + 69*a0)*\[Omega]^3) - 90*M*(4*a0^4*\[Omega]^2 - 
            4*a0^3*\[Omega]*(3*I + \[Omega]) + a0^2*(-12 + (15*I)*\[Omega] - 
              8*\[Omega]^2) - 4*(2 - (3*I)*\[Omega] + 8*\[Omega]^2) + 
            a0*(-1 - (12*I)*\[Omega] + 16*\[Omega]^2)) + 
          15*M^2*(16*a0^3*(12 + I*\[Omega])*\[Omega]^2 + 
            (6*I)*a0^4*\[Omega]^3 + a0^2*\[Omega]*(-502*I + 129*\[Omega]) - 
            2*a0*(93 + (56*I)*\[Omega] + 190*\[Omega]^2) + 
            4*(-12 + (14*I)*\[Omega] + 65*\[Omega]^2 + (8*I)*\[Omega]^3)))))/
      (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (I*ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*((-6*I)*(4 + n + (2*I)*\[Omega]) + 
            M^2*\[Omega]*(12 - (8*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^
                2) + M*(-24*I - (6*I)*n - 2*(-10 + a0)*\[Omega] + 
              (2 + a0)*n*\[Omega] + (2*I)*(-4 + 4*a0 + 3*a0^2)*\[Omega]^2)) - 
          ang*(-1 + \[ScriptM])*\[Omega]*
           (n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
            \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - 
                (4*I)*\[Omega] + (3*I)*a0^2*\[Omega]) + M^2*(12 - 
                (8*I)*(2 + a0)*\[Omega] + (2 + a0)^2*\[Omega]^2))))*
         SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] - 
        ang*\[Omega]*(n*(-6*I + M*(-6*I + (2 + a0)*\[Omega])) + 
          \[Omega]*(12 + 2*M*(6 + a0*(-3 + (4*I)*\[Omega]) - (4*I)*\[Omega] + 
              (3*I)*a0^2*\[Omega]) + M^2*(12 - (8*I)*(2 + a0)*\[Omega] + 
              (2 + a0)^2*\[Omega]^2)))*Conjugate[
          1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (\[Omega]*(-((2 + a0)*M*\[Omega]*
          ((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
                I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*
                 \[Omega] - I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 
                5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
                 \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
                  9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
                 (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
                 \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                  77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                  (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                   \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                  (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                   \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                  2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
                5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 
                    90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 
                    8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 
                    2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
           (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                 (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 
                    3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 
                    3*M)*M*\[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                       \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                    (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                  \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                      3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                      2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*
               SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
              ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
                I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*
                 (-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
                 (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
                1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*Sqrt[
                (2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + \[ScriptL]*
                   (3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
               SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
            (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
             \[Omega]^2))) + (3*I)*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3) + 2*(1 + M)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*
                (-I + 8*(1 + M)*\[Omega] + I*(2 + a0)*M*\[Omega]^2) + 10*n^2*
                (5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 24*M^2)*
                  \[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 5*n*
                (12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
                  \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
                   9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
                  (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
                  \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                   77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                   (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                    \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                   (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                    \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                   2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
                 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 
                     90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 
                     8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 
                     2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
            (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                  (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 
                     3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 
                     3*M)*M*\[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                        \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                     (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                   \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                       3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                       2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*
                SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + ang*
                \[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
                 I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 2*(2 + a0)*
                  (-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
                  (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
                 1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
                Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                   \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
                SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
             (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
              \[Omega]^2)))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/((-2 + r)*r^2*
    Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]), 
 I*K\[Infinity]0*r*\[Omega] - K\[Infinity]0*(n - (2*I)*(1 + M)*\[Omega]) + 
  ((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
         I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
         I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
         (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 
       5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
          \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
           9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - 
         I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
          \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
           77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*
          (-360 - (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
            \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
           (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
            \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
           2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
         5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
            \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*\[Omega]^3 + 
           12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^4))))/
     (120*\[Omega]^3) + (ang*Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(36*I + I*n^2 + 
           81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 
           2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
           n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
         2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
             (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
             2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
             M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
         Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
         I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 
         2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
         n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^2))/
   r^3 + 
  ((K\[Infinity]0*((-3*I)*n^2 - (3*I)*n*(1 - (4*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2) + \[Omega]*(-9 + (24*I)*\[Omega] - 
         2*M^2*\[Omega]*(-12*I + 5*(2 + a0)*\[Omega]) + 
         3*M*(-3 - (4*I)*(-2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
            \[Omega]^2))))/(6*\[Omega]) - 
    (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
          \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) - 
    (I/2)*\[Omega]*
     (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)) + 
      (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
           21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
              \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
              \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
             (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
             2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
           15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
             2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
             2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
             a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^
                 3)))))/(30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
        ((-I)*\[ScriptM]*(2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - 
             (9*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*
              \[Omega]^3) + ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
             54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
           \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
          (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
           (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
          Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
          Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
             \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
          SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
       (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
        (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
      (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
             18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
            \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*
                M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                   2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
           0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
              \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
               \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
             \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
            (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
           2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
        n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2))))/r + 
  ((K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
         21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
          \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
         27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
            \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
         3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
            \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
       I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
         2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*\[Omega] + 
           (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
         5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
           (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
           2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
         15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
           2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
           2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
           a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
     (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
        (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
            \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
         ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
           54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
           2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
         \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
        (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
         (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
        Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
      (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 
            5*n^3*(-I + 8*(1 + M)*\[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
            10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^
                3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + 
                (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + (142 - 85*a0 + 
                  66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - 
              I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*
               \[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 
                77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*(-360 - 
                (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 653*a0^2)*
                 \[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - 
                (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*
                 \[Omega] - (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
                2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
              5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
                 \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*
                 \[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*
                 \[Omega]^4))))/(120*\[Omega]^3) + 
         (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(
                36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                  16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
              2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
                  (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
                  2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
                  M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 
              81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^
                2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
              n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
             Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
             Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
             SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
          (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
           \[Omega]^2)) + \[Omega]*
        ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
            10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
              I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
            5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
              (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
              M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                 \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
            I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
              2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
              5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                  19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                 \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
          (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
           (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                   \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*((3*I)*n^2 + 
                n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
                \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                  M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + 
              n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
              \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                   \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                  \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
                \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - 
                 \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], 
              Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*
           \[Sigma]*\[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
      (2 + a0)*M*\[Omega]^2))/r^2, I*K\[Infinity]0*\[Omega] - 
  (E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*K\[Infinity]0*
    r^2*(a0 + r)*\[Omega]^2)/((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + 
      r^2]) - 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*K\[Infinity]0*
    r*(a0 + r)*\[Omega]*(n - (2*I)*(1 + M)*\[Omega]))/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (3*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
          I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
          I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
          (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 
        5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
           \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
            9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
           (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
        \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 
          2*M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
            (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 
          15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 
          5*M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - 
            (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
            2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
          5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
             \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*\[Omega]^3 + 
            12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^4))))/
      (120*\[Omega]^3) + (ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(36*I + I*n^2 + 
            81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*
             \[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
            n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
          2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
              (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
              2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
              M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
          I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 
          2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
          n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
       \[Omega]^2)))/r^4 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
          I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
          I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
          (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 
        5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
           \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 
            9*(-6 + 29*a0)*M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*
           (24 + 18*a0^2 + 18*M + 52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
        \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 
          2*M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
            (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 
          15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 
          5*M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - 
            (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 
            2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) - 
          5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
             \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*\[Omega]^3 + 
            12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^4))))/
      (120*\[Omega]^3) + (ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(36*I + I*n^2 + 
            81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*
             \[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
            n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
          2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
              (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
              2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
              M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
          I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 
          2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
          n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
       \[Omega]^2)))/((-2 + r)*r^2*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + 
      r^2]) - 
  ((K\[Infinity]0*((-3*I)*n^2 - (3*I)*n*(1 - (4*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2) + \[Omega]*(-9 + (24*I)*\[Omega] - 
         2*M^2*\[Omega]*(-12*I + 5*(2 + a0)*\[Omega]) + 
         3*M*(-3 - (4*I)*(-2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
            \[Omega]^2))))/(6*\[Omega]) - 
    (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
          \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) - 
    (I/2)*\[Omega]*
     (-1/30*(K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)) + 
      (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
           21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
              \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
              \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
             (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
             2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
           15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
             2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
             2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
             a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^
                 3)))))/(30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
        ((-I)*\[ScriptM]*(2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - 
             (9*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*
              \[Omega]^3) + ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
             54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
           \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
          (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
           (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
          Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
          Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
             \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
          SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
       (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
        (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
      (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
             18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
            \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*
                M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                   2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
           0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
              \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
               \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
             \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
            (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
           2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
        n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2))))/r^2 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*((-3*I)*n^2 - 
        (3*I)*n*(1 - (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
        \[Omega]*(-9 + (24*I)*\[Omega] - 2*M^2*\[Omega]*
           (-12*I + 5*(2 + a0)*\[Omega]) + 
          3*M*(-3 - (4*I)*(-2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
             \[Omega]^2))))/(6*\[Omega]) - 
     (I*ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
           \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
        ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]) - 
     (I/2)*\[Omega]*(-1/30*(K\[Infinity]0*(45*n^3 - 
           5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 21*(2 + a0)*M*\[Omega]^2 + 
             (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
           5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 27*(3*I + 8*\[Omega]) + 
             M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*\[Omega] + (20 - 92*a0 - 
                 51*a0^2)*\[Omega]^2) + 3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - 
               (2*I)*(-12 + 8*a0 + 7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 
                 3*a0^2)*\[Omega]^3)) + I*\[Omega]^2*
            (-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
             2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*\[Omega] + 
               (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
             5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                 9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
                \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
               8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*
                a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + a0*
                (-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
         (\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
            21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
               \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 
                7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
               \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - (9*I)*(-36 + 28*a0 + 
                23*a0^2)*\[Omega]^2 + 2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*
               \[Omega]^3) + 15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(
                -3*I + 2*\[Omega]) + 2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 
                4*\[Omega]^2) + 2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
              a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
           (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
               \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
            ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
              54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
            \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
           (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
            (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
           Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
           Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
              \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
           SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
        (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
       (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
              18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
             \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
                3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                    2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
            0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
               \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                 \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
              \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
             (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
            2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*
         enr*n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)))))/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (2*((K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
          21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
           \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
          27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
             \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
          3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
             \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
        I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
          2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*\[Omega] + 
            (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
          5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
            (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
            2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
          15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
            2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
            2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
            a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
      (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
         (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
             \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
          ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
            54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
            2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
          \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
         (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
          (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*
                \[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
             10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*
                \[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*
                (8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + 
                 (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*
                \[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - 
                 a0*(24 + 71*M))*\[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*
                (-240*I + 77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*
                (-360 - (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 
                   653*a0^2)*\[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 
                 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*
                (33*I + (22 - 205*a0)*\[Omega] - (4*I)*(54 - 53*a0 + 
                   104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) - 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - 
                 (142 + 35*a0 + 90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 
                   19*a0^2 + 8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + 
                   a0^3 + 2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
          (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                   16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                  \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                      \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                   (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                 \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                     3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                     2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*SphericalHarmonicY[
               \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
              (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^2)) + \[Omega]*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/r^3 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
          21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
           \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
          27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
             \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
          3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
             \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
        I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 270*(3*I + 8*\[Omega]) + 
          2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*\[Omega] + 
            (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
          5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
            (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
            2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
          15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
            2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
            2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
            a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
      (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
         (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
             \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
          ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
            54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
            2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
          \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
         (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
          (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
         Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (\[Omega]*((-3*I)*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*
                \[Omega] + I*(2 + a0)*M*\[Omega]^2) + 
             10*n^2*(5*I + 11*(1 + M)*\[Omega] - I*(24 + 14*M - 17*a0*M + 
                 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 3*a0 - 5*M)*M*
                \[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*
                (8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - (192 + 
                 (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 192*M^3)*
                \[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 52*M^2 - 
                 a0*(24 + 71*M))*\[Omega]^4) - \[Omega]*(8*M^4*\[Omega]^3*
                (-240*I + 77*(2 + a0)*\[Omega]) - 2*M^3*\[Omega]^2*
                (-360 - (5*I)*(-2 + 383*a0)*\[Omega] + (12 + 1312*a0 + 
                   653*a0^2)*\[Omega]^2) + 15*(-12 + (11*I)*\[Omega] + 
                 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 5*M^2*\[Omega]*
                (33*I + (22 - 205*a0)*\[Omega] - (4*I)*(54 - 53*a0 + 
                   104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) - 5*M*(36 - (6*I)*(19 + 4*a0)*\[Omega] - 
                 (142 + 35*a0 + 90*a0^2)*\[Omega]^2 - (6*I)*(-58 + 29*a0 - 
                   19*a0^2 + 8*a0^3)*\[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + 
                   a0^3 + 2*a0^4)*\[Omega]^4))))/(120*\[Omega]^3) + 
          (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*
                (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                   16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                  \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*
                      \[Omega]))) + 2*enr*n*((-I)*n^2 - I*n*(2 - 
                   (4*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2) + 
                 \[Omega]*(-49 + (16*I)*\[Omega] - 2*M^2*\[Omega]*(-8*I + 
                     3*(2 + a0)*\[Omega]) + M*(-49 + I*(34 + a0)*\[Omega] + 
                     2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2))))*SphericalHarmonicY[
               \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
              (36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                 16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^2)) + \[Omega]*
         ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
             10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + I*
                (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 5*n*\[Omega]*
              (9*I + 52*(2 + a0)*M^3*\[Omega]^3 - (2 + a0)*M^2*\[Omega]^2*
                (11*I + (-18 + 71*a0)*\[Omega]) + M*(9*I + 5*(2 + a0)*
                  \[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*\[Omega]^2 + 
                 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
             I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 2*(2 + a0)*
                M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 5*M^2*
                (-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
                  \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                  \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                 (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                 4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
           (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
            (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                 (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*
                  M*\[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                    \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
                ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*
                      \[Omega])) + \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                   M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                      \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
               \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
                (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
                (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                    \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
              Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
              Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                 \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
              SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
           (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
            \[Omega]^3))))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/((-2 + r)*r*
    Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]), 
 K\[Infinity]0 + (K\[Infinity]0*(-2 - 2*M - (I*n)/\[Omega]) + 
    (I*K\[Infinity]0*(n - (2*I)*(1 + M)*\[Omega]))/\[Omega])/r + 
  ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
       10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
         I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
       5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
         (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
         M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
            \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
       I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
         2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
         5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
            \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) + 
         15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
           (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
           4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
     (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
           (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
            \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
              \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
          ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
           \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
             M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
         Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
          (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
         \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
           M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
              \[Omega]^2)))*Conjugate[
         1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^3))/
   r^4 + ((K\[Infinity]0*(3*n^2 - 3*n*(-1 + (2 + a0)*M*\[Omega]^2) - 
       I*\[Omega]*(9 - 10*(2 + a0)*M^2*\[Omega]^2 + 
         3*M*(3 + (2*I)*(2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*\[Omega]^2))))/
     (6*\[Omega]^2) - (ang*Sqrt[2*Pi]*\[Mu]*
      (\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*\[Omega])*
        SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
       ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]) + 
    ((K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
           21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
              \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
              \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - 
             (15*I)*(4 + 20*a0 + 9*a0^2)*\[Omega]^2 + 
             2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*\[Omega]^3) + 
           15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 8*a0^4*\[Omega]^3 + 
             2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*
              (-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
             a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 32*\[Omega]^3)))))/
       (30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
      (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
           21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
            \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
           27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
              \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
           3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 7*a0^2)*
              \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
         I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
           270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
              \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
           5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - 
             (9*I)*(-36 + 28*a0 + 23*a0^2)*\[Omega]^2 + 
             2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*\[Omega]^3) + 
           15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) + 
             2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 4*\[Omega]^2) + 
             2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
             a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 32*\[Omega]^
                 3)))))/(30*\[Omega]^2*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)) - (ang*Sqrt[Pi/2]*\[Mu]*
        ((-I)*\[ScriptM]*(2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - 
             (9*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*
              \[Omega]^3) + ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
             54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
           \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
          (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
           (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
          Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
          Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
             \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
          SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
       (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
        (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
      (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
             18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
             2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
            \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 3*
                M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                   2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
           0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
              \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
               \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
             \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
            (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
           2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*
        n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
         (2 + a0)*M*\[Omega]^2)))/2)/r^2 + 
  ((K\[Infinity]0*(15*n^3*(6*I + 6*(1 + M)*\[Omega] + 
         I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(9*I + 9*(1 + M)*\[Omega] - 
         (6*I)*(3 + 2*M - 2*a0*M + 3*M^2)*\[Omega]^2 - 
         (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - (5*I)*n*\[Omega]*
        (4*M^3*\[Omega]^2*(-36*I + 13*(2 + a0)*\[Omega]) + 
         18*(3*I + 5*\[Omega] - (8*I)*\[Omega]^2) + M^2*\[Omega]*
          (90 + I*(14 + 223*a0)*\[Omega] - 5*(-4 + 28*a0 + 15*a0^2)*
            \[Omega]^2) + 3*M*(18*I + (6 - 27*a0)*\[Omega] - 
           I*(14 - 37*a0 + 14*a0^2)*\[Omega]^2 + 2*(8 - 4*a0 + 2*a0^2 + 
             3*a0^3)*\[Omega]^3)) + \[Omega]^2*
        (180*(-3 + (8*I)*\[Omega])*\[Omega] - 8*M^4*\[Omega]^2*
          (-180*I + 77*(2 + a0)*\[Omega]) + 2*M^3*\[Omega]*
          (-270 - (15*I)*(14 + 103*a0)*\[Omega] + (92 + 1392*a0 + 673*a0^2)*
            \[Omega]^2) - 5*M^2*\[Omega]*(-234 - (84*I)*\[Omega] + 
           88*\[Omega]^2 + 168*a0^3*\[Omega]^2 + 3*a0^2*\[Omega]*
            (-125*I + 98*\[Omega]) + a0*(-279 + (72*I)*\[Omega] - 
             40*\[Omega]^2)) + 15*M*(48*I - 2*\[Omega] + (84*I)*\[Omega]^2 - 
           32*\[Omega]^3 + 8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*
            (-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*(30 - (33*I)*\[Omega] + 
             4*\[Omega]^2) + a0*(24*I - 67*\[Omega] + (6*I)*\[Omega]^2 + 
             16*\[Omega]^3)))))/(30*\[Omega]^3*(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
      (\[ScriptM]*(2*enr*n*(-72*I - 144*(1 + M)*\[Omega] + 
           (3*I)*(12 + (26 + a0)*M + 12*M^2)*\[Omega]^2 + 
           2*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3 - 
           (3*I)*n*(6 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
         ang*(-1 + \[ScriptM])*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
             I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
             2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
             3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
                \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
         Pi/2, 0] + ang*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
           I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
           2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
           3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
              \[Omega]^2)))*Conjugate[
         1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
        Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
           \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
        SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
     (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^2*
      (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
    (-3*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
             I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
             I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - 
             (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) + 
           5*n*(12*I + 5*(1 + M)*\[Omega] - (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*
              \[Omega]^2 - (192 + (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*
                M^2 + 192*M^3)*\[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*
                M + 52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
           \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 
             2*M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
               (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 
             15*(-12 + (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*
                \[Omega]^3) + 5*M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - 
               (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 
                 123*a0^2 + 78*a0^3)*\[Omega]^3) - 5*M*(36 - (6*I)*
                (19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*\[Omega]^2 - 
               (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*\[Omega]^3 + 12*
                (-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^4))))/
         (120*\[Omega]^3) + (ang*Sqrt[Pi/2]*\[Mu]*
          (\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(36*I + I*n^2 + 81*
                (1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*
                \[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*
                \[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
             2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
                 (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
                 2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
                 M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                    \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
             Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 81*(1 + M)*\[Omega] - 
             I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^2 - 
             2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + n*\[Omega]*
              (4 + M*(4 + I*(2 + a0)*\[Omega])))*Conjugate[
             1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
            Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + \[ScriptL]*
                (3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
            SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
         (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
          \[Omega]^2)) - I*\[Omega]*
       ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
           10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
             I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
           5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
             (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
             M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
           I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
             2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
             5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                 19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
               (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 4*(-8 + 4*a0 - 
                 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/(120*\[Omega]^4) + 
        ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*
            (2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - (37*I)*
                (2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                  \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
              ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
               \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                 M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                    \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
             Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
              (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + \[Omega]*
              (27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + (39*I)*(2 + a0)*
                  \[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*\[Omega]^2)))*
            Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
            Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + \[ScriptL]*
                (3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
            SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
          \[Omega]^3)))/(9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2))/
   r^3, 
 (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*K\[Infinity]0*
    r*(a0 + r)*\[Omega])/((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + 
      r^2]) - (K\[Infinity]0*(-2 - 2*M - (I*n)/\[Omega]) + 
    (I*K\[Infinity]0*(n - (2*I)*(1 + M)*\[Omega]))/\[Omega])/r^2 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*(K\[Infinity]0*(-2 - 2*M - (I*n)/\[Omega]) + 
     (I*K\[Infinity]0*(n - (2*I)*(1 + M)*\[Omega]))/\[Omega]))/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (4*((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
        10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
          I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
        5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
          (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
          M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
             \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
        I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
          2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
          5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
             \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) + 
          15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
            (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
            4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
      (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
            (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
             \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
               \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
           ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
            \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + 
                (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
           (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
          \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
            M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
               \[Omega]^2)))*Conjugate[
          1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^3)))/
   r^5 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
        10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
          I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
        5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
          (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
          M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
             \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
        I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
          2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
          5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 19*a0^2)*
             \[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*\[Omega]^3) + 
          15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
            (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
            4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
      (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
            (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
             \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
               \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*
           ((3*I)*n^2 + n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
            \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + M*(27 + 
                (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + n*\[Omega]*
           (-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
          \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
            M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
               \[Omega]^2)))*Conjugate[
          1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^3)))/
   ((-2 + r)*r^3*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (2*((K\[Infinity]0*(3*n^2 - 3*n*(-1 + (2 + a0)*M*\[Omega]^2) - 
        I*\[Omega]*(9 - 10*(2 + a0)*M^2*\[Omega]^2 + 
          3*M*(3 + (2*I)*(2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
             \[Omega]^2))))/(6*\[Omega]^2) - 
     (ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
           \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
        ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]) + 
     ((K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
            21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
               \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 
                7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
               \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - (9*I)*(-36 + 28*a0 + 
                23*a0^2)*\[Omega]^2 + 2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*
               \[Omega]^3) + 15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(
                -3*I + 2*\[Omega]) + 2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 
                4*\[Omega]^2) + 2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
              a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
           (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
               \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
            ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
              54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
            \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
           (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
            (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
           Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
           Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
              \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
           SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
        (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
              18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
             \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
                3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                    2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
            0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
               \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                 \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
              \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
             (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
            2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*
         enr*n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)))/2))/r^3 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*(3*n^2 - 3*n*(-1 + (2 + a0)*M*\[Omega]^2) - 
        I*\[Omega]*(9 - 10*(2 + a0)*M^2*\[Omega]^2 + 
          3*M*(3 + (2*I)*(2 + a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
             \[Omega]^2))))/(6*\[Omega]^2) - 
     (ang*Sqrt[2*Pi]*\[Mu]*(\[ScriptM]*(-2*enr*n + ang*(-1 + \[ScriptM])*
           \[Omega])*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 0] + 
        ang*\[Omega]*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]) + 
     ((K\[Infinity]0*(45*n^3 - 5*n^2*(-9 + (54*I)*(1 + M)*\[Omega] + 
            21*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) + M^2*\[Omega]*(-216 + (64*I)*(2 + a0)*
               \[Omega] + (20 - 92*a0 - 51*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-6 + 5*a0)*\[Omega] - (2*I)*(-12 + 8*a0 + 7*a0^2)*
               \[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (240*I)*(2 + a0)*
               \[Omega] + (-28 + 912*a0 + 463*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 27*(-2 + 23*a0)*\[Omega] - (15*I)*(4 + 20*a0 + 
                9*a0^2)*\[Omega]^2 + 2*(20 - 44*a0 + 105*a0^2 + 66*a0^3)*
               \[Omega]^3) + 15*M*(12*I + 82*\[Omega] - (48*I)*\[Omega]^2 + 
              8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(-3*I + 2*\[Omega]) - 
              2*a0^2*\[Omega]*(-24 + (3*I)*\[Omega] + 4*\[Omega]^2) + 
              a0*(-48*I - 79*\[Omega] - (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (K\[Infinity]0*(45*n^3 + 5*n^2*(9 - (54*I)*(1 + M)*\[Omega] - 
            21*(2 + a0)*M*\[Omega]^2 - (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
             \[Omega]^3) + 5*n*\[Omega]*(28*(2 + a0)*M^3*\[Omega]^3 - 
            27*(3*I + 8*\[Omega]) - M^2*\[Omega]*(216 - (100*I)*(2 + a0)*
               \[Omega] + (4 + 116*a0 + 57*a0^2)*\[Omega]^2) + 
            3*M*(-27*I + 9*(-10 + 3*a0)*\[Omega] - (2*I)*(-24 + 2*a0 + 
                7*a0^2)*\[Omega]^2 + 2*a0*(-8 + 2*a0 + 3*a0^2)*\[Omega]^3)) + 
          I*\[Omega]^2*(-376*(2 + a0)*M^4*\[Omega]^3 + 
            270*(3*I + 8*\[Omega]) + 2*M^3*\[Omega]*(1080 - (540*I)*(2 + a0)*
               \[Omega] + (172 + 1112*a0 + 513*a0^2)*\[Omega]^2) - 
            5*M^2*(-162*I + 9*(-46 + 49*a0)*\[Omega] - (9*I)*(-36 + 28*a0 + 
                23*a0^2)*\[Omega]^2 + 2*(-4 - 44*a0 + 123*a0^2 + 72*a0^3)*
               \[Omega]^3) + 15*M*(8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*(
                -3*I + 2*\[Omega]) + 2*a0^2*\[Omega]*(6 + (9*I)*\[Omega] - 
                4*\[Omega]^2) + 2*(6*I + 77*\[Omega] - (48*I)*\[Omega]^2) + 
              a0*(-48*I - 115*\[Omega] + (12*I)*\[Omega]^2 + 
                32*\[Omega]^3)))))/(30*\[Omega]^2*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (ang*Sqrt[Pi/2]*\[Mu]*((-I)*\[ScriptM]*
           (2*enr*n*(18*I - (9*I)*n + 54*(1 + M)*\[Omega] - (9*I)*(2 + a0)*M*
               \[Omega]^2 - 2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + 
            ang*(-1 + \[ScriptM])*\[Omega]*(-54*I + (9*I)*n - 
              54*(1 + M)*\[Omega] + (9*I)*(2 + a0)*M*\[Omega]^2 + 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3))*SphericalHarmonicY[
            \[ScriptL], \[ScriptM], Pi/2, 0] + ang*\[Omega]*
           (-54 + 9*n + (54*I)*(1 + M)*\[Omega] + 9*(2 + a0)*M*\[Omega]^2 - 
            (2*I)*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3)*
           Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
           Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
              \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
           SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
        (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]*
         (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) - 
       (ang*Sqrt[Pi/2]*\[Mu]*(I*\[ScriptM]*(2*enr*n*(-36*I - (9*I)*n + 
              18*(1 + M)*\[Omega] - (15*I)*(2 + a0)*M*\[Omega]^2 - 
              2*(2 + a0)*(3*a0 - 2*M)*M*\[Omega]^3) + ang*(-1 + \[ScriptM])*
             \[Omega]*((9*I)*n + \[Omega]*(-18 - 4*(2 + a0)*M^2*\[Omega]^2 + 
                3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*\[Omega]^
                    2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], Pi/2, 
            0] - ang*\[Omega]*(9*n - I*\[Omega]*(-18 - 4*(2 + a0)*M^2*
               \[Omega]^2 + 3*M*(-6 + (5*I)*(2 + a0)*\[Omega] + 2*a0*(2 + a0)*
                 \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
              \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*
             (\[ScriptL] - \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 
            2 + \[ScriptM], Pi/2, 0]))/(3*E^((r - rpart)^2/(2*\[Sigma]^2))*
         enr*n*(1 + n)*\[Sigma]*\[Omega]*(9 - (6*I)*(1 + M)*\[Omega] + 
          (2 + a0)*M*\[Omega]^2)))/2))/
   ((-2 + r)*r*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (3*((K\[Infinity]0*(15*n^3*(6*I + 6*(1 + M)*\[Omega] + 
          I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(9*I + 9*(1 + M)*\[Omega] - 
          (6*I)*(3 + 2*M - 2*a0*M + 3*M^2)*\[Omega]^2 - 
          (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - (5*I)*n*\[Omega]*
         (4*M^3*\[Omega]^2*(-36*I + 13*(2 + a0)*\[Omega]) + 
          18*(3*I + 5*\[Omega] - (8*I)*\[Omega]^2) + M^2*\[Omega]*
           (90 + I*(14 + 223*a0)*\[Omega] - 5*(-4 + 28*a0 + 15*a0^2)*
             \[Omega]^2) + 3*M*(18*I + (6 - 27*a0)*\[Omega] - 
            I*(14 - 37*a0 + 14*a0^2)*\[Omega]^2 + 2*(8 - 4*a0 + 2*a0^2 + 
              3*a0^3)*\[Omega]^3)) + \[Omega]^2*
         (180*(-3 + (8*I)*\[Omega])*\[Omega] - 8*M^4*\[Omega]^2*
           (-180*I + 77*(2 + a0)*\[Omega]) + 2*M^3*\[Omega]*
           (-270 - (15*I)*(14 + 103*a0)*\[Omega] + (92 + 1392*a0 + 673*a0^2)*
             \[Omega]^2) - 5*M^2*\[Omega]*(-234 - (84*I)*\[Omega] + 
            88*\[Omega]^2 + 168*a0^3*\[Omega]^2 + 3*a0^2*\[Omega]*
             (-125*I + 98*\[Omega]) + a0*(-279 + (72*I)*\[Omega] - 
              40*\[Omega]^2)) + 15*M*(48*I - 2*\[Omega] + (84*I)*\[Omega]^2 - 
            32*\[Omega]^3 + 8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*
             (-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*(30 - (33*I)*\[Omega] + 
              4*\[Omega]^2) + a0*(24*I - 67*\[Omega] + (6*I)*\[Omega]^2 + 
              16*\[Omega]^3)))))/(30*\[Omega]^3*(9 - (6*I)*(1 + M)*\[Omega] + 
        (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*(-72*I - 144*(1 + M)*\[Omega] + 
            (3*I)*(12 + (26 + a0)*M + 12*M^2)*\[Omega]^2 + 
            2*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3 - 
            (3*I)*n*(6 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
          ang*(-1 + \[ScriptM])*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
              I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
              2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
              3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
            I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
            2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
            3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*\[Omega]^
                2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^2*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (-3*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
              I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
              I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 
                3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - 
              (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - 
              (192 + (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 
                192*M^3)*\[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 
                52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
            \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 
              2*M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
                (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 15*(-12 + 
                (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 
              5*M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - 
                (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 
                  123*a0^2 + 78*a0^3)*\[Omega]^3) - 5*M*(36 - 
                (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
                 \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*
                 \[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*
                 \[Omega]^4))))/(120*\[Omega]^3) + 
         (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(
                36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                  16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
              2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
                  (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
                  2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
                  M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 
              81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^
                2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
              n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
             Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
             Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
             SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
          (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
           \[Omega]^2)) - I*\[Omega]*
        ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
            10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
              I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
            5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
              (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
              M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                 \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
            I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
              2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
              5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                  19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                 \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
          (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
           (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                   \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*((3*I)*n^2 + 
                n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
                \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                  M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + 
              n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
              \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                   \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                  \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
                \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - 
                 \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], 
              Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*
           \[Sigma]*\[Omega]^3)))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/r^4 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*
       (Pi - 2*ArcTan[(a0 - M + r)/Sqrt[(4 + 2*a0 - M)*M]]))/2)*(a0 + r)*
    \[Omega]*((K\[Infinity]0*(15*n^3*(6*I + 6*(1 + M)*\[Omega] + 
          I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(9*I + 9*(1 + M)*\[Omega] - 
          (6*I)*(3 + 2*M - 2*a0*M + 3*M^2)*\[Omega]^2 - 
          (2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - (5*I)*n*\[Omega]*
         (4*M^3*\[Omega]^2*(-36*I + 13*(2 + a0)*\[Omega]) + 
          18*(3*I + 5*\[Omega] - (8*I)*\[Omega]^2) + M^2*\[Omega]*
           (90 + I*(14 + 223*a0)*\[Omega] - 5*(-4 + 28*a0 + 15*a0^2)*
             \[Omega]^2) + 3*M*(18*I + (6 - 27*a0)*\[Omega] - 
            I*(14 - 37*a0 + 14*a0^2)*\[Omega]^2 + 2*(8 - 4*a0 + 2*a0^2 + 
              3*a0^3)*\[Omega]^3)) + \[Omega]^2*
         (180*(-3 + (8*I)*\[Omega])*\[Omega] - 8*M^4*\[Omega]^2*
           (-180*I + 77*(2 + a0)*\[Omega]) + 2*M^3*\[Omega]*
           (-270 - (15*I)*(14 + 103*a0)*\[Omega] + (92 + 1392*a0 + 673*a0^2)*
             \[Omega]^2) - 5*M^2*\[Omega]*(-234 - (84*I)*\[Omega] + 
            88*\[Omega]^2 + 168*a0^3*\[Omega]^2 + 3*a0^2*\[Omega]*
             (-125*I + 98*\[Omega]) + a0*(-279 + (72*I)*\[Omega] - 
              40*\[Omega]^2)) + 15*M*(48*I - 2*\[Omega] + (84*I)*\[Omega]^2 - 
            32*\[Omega]^3 + 8*a0^4*\[Omega]^3 + 2*a0^3*\[Omega]^2*
             (-3*I + 2*\[Omega]) - 2*a0^2*\[Omega]*(30 - (33*I)*\[Omega] + 
              4*\[Omega]^2) + a0*(24*I - 67*\[Omega] + (6*I)*\[Omega]^2 + 
              16*\[Omega]^3)))))/(30*\[Omega]^3*(9 - (6*I)*(1 + M)*\[Omega] + 
        (2 + a0)*M*\[Omega]^2)) + (ang*Sqrt[Pi/2]*\[Mu]*
       (\[ScriptM]*(2*enr*n*(-72*I - 144*(1 + M)*\[Omega] + 
            (3*I)*(12 + (26 + a0)*M + 12*M^2)*\[Omega]^2 + 
            2*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3 - 
            (3*I)*n*(6 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
          ang*(-1 + \[ScriptM])*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
              I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
              2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
              3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*
                 \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], \[ScriptM], 
          Pi/2, 0] + ang*\[Omega]*(3*n*(6*I + 6*(1 + M)*\[Omega] + 
            I*(2 + a0)*M*\[Omega]^2) + \[Omega]*(72 - (36*I)*\[Omega] + 
            2*M^2*\[Omega]*(-18*I + 5*(2 + a0)*\[Omega]) - 
            3*M*(-24 + I*(34 + 5*a0)*\[Omega] + 2*(-2 + a0 + a0^2)*\[Omega]^
                2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
         Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
            \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
         SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
      (3*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*\[Omega]^2*
       (9 - (6*I)*(1 + M)*\[Omega] + (2 + a0)*M*\[Omega]^2)) + 
     (-3*((K\[Infinity]0*((5*I)*n^4 + 5*n^3*(-I + 8*(1 + M)*\[Omega] + 
              I*(2 + a0)*M*\[Omega]^2) + 10*n^2*(5*I + 11*(1 + M)*\[Omega] - 
              I*(24 + 14*M - 17*a0*M + 24*M^2)*\[Omega]^2 - (2 + a0)*(-3 + 
                3*a0 - 5*M)*M*\[Omega]^3) + 5*n*(12*I + 5*(1 + M)*\[Omega] - 
              (9*I)*(8 + (2 - 7*a0)*M + 8*M^2)*\[Omega]^2 - 
              (192 + (142 - 85*a0 + 66*a0^2)*M - 9*(-6 + 29*a0)*M^2 + 
                192*M^3)*\[Omega]^3 - I*(2 + a0)*M*(24 + 18*a0^2 + 18*M + 
                52*M^2 - a0*(24 + 71*M))*\[Omega]^4) - 
            \[Omega]*(8*M^4*\[Omega]^3*(-240*I + 77*(2 + a0)*\[Omega]) - 
              2*M^3*\[Omega]^2*(-360 - (5*I)*(-2 + 383*a0)*\[Omega] + 
                (12 + 1312*a0 + 653*a0^2)*\[Omega]^2) + 15*(-12 + 
                (11*I)*\[Omega] + 48*\[Omega]^2 - (128*I)*\[Omega]^3) + 
              5*M^2*\[Omega]*(33*I + (22 - 205*a0)*\[Omega] - 
                (4*I)*(54 - 53*a0 + 104*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 
                  123*a0^2 + 78*a0^3)*\[Omega]^3) - 5*M*(36 - 
                (6*I)*(19 + 4*a0)*\[Omega] - (142 + 35*a0 + 90*a0^2)*
                 \[Omega]^2 - (6*I)*(-58 + 29*a0 - 19*a0^2 + 8*a0^3)*
                 \[Omega]^3 + 12*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*
                 \[Omega]^4))))/(120*\[Omega]^3) + 
         (ang*Sqrt[Pi/2]*\[Mu]*(\[ScriptM]*(ang*(-1 + \[ScriptM])*\[Omega]*(
                36*I + I*n^2 + 81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 
                  16*M^2)*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega]))) + 
              2*enr*n*((-I)*n^2 - I*n*(2 - (4*I)*(1 + M)*\[Omega] + 
                  (2 + a0)*M*\[Omega]^2) + \[Omega]*(-49 + (16*I)*\[Omega] - 
                  2*M^2*\[Omega]*(-8*I + 3*(2 + a0)*\[Omega]) + 
                  M*(-49 + I*(34 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*(36*I + I*n^2 + 
              81*(1 + M)*\[Omega] - I*(16 + 38*M + 3*a0*M + 16*M^2)*\[Omega]^
                2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*\[Omega]^3 + 
              n*\[Omega]*(4 + M*(4 + I*(2 + a0)*\[Omega])))*
             Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - \[ScriptM]]]]*
             Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + \[ScriptM]^2 + 
                \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - \[ScriptM])!]*
             SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], Pi/2, 0]))/
          (12*E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*\[Sigma]*
           \[Omega]^2)) - I*\[Omega]*
        ((K\[Infinity]0*(-15*n^4 + 5*n^3*(3 + (2 + a0)*M*\[Omega]^2) + 
            10*n^2*(3 + (9*I)*(1 + M)*\[Omega] - 5*(2 + a0)*M*\[Omega]^2 + 
              I*(2 + a0)*(-3 + 3*a0 - 5*M)*M*\[Omega]^3) - 
            5*n*\[Omega]*(9*I + 52*(2 + a0)*M^3*\[Omega]^3 - 
              (2 + a0)*M^2*\[Omega]^2*(11*I + (-18 + 71*a0)*\[Omega]) + 
              M*(9*I + 5*(2 + a0)*\[Omega] + (3*I)*(6 + 23*a0 + 10*a0^2)*
                 \[Omega]^2 + 6*(8 - 4*a0 + 2*a0^2 + 3*a0^3)*\[Omega]^3)) - 
            I*\[Omega]^2*(135*I - 616*(2 + a0)*M^4*\[Omega]^3 + 
              2*(2 + a0)*M^3*\[Omega]^2*(45*I + (6 + 653*a0)*\[Omega]) - 
              5*M^2*(-27*I + 15*(2 + a0)*\[Omega] + (4*I)*(22 + 49*a0 + 
                  19*a0^2)*\[Omega]^2 + 2*(44 - 44*a0 + 123*a0^2 + 78*a0^3)*
                 \[Omega]^3) + 15*M*(18*I - 3*(6 + 7*a0 + 2*a0^2)*\[Omega] + 
                (6*I)*(-2 + a0 + 9*a0^2 + 4*a0^3)*\[Omega]^2 + 
                4*(-8 + 4*a0 - 2*a0^2 + a0^3 + 2*a0^4)*\[Omega]^3))))/
          (120*\[Omega]^4) + ((I/12)*ang*Sqrt[Pi/2]*\[Mu]*
           (\[ScriptM]*(2*enr*n*(36*I - (3*I)*n^2 + 117*(1 + M)*\[Omega] - 
                (37*I)*(2 + a0)*M*\[Omega]^2 - 2*(2 + a0)*(-1 + 3*a0 - 3*M)*M*
                 \[Omega]^3 + n*(-6*I + 12*(1 + M)*\[Omega] + I*(2 + a0)*M*
                   \[Omega]^2)) + ang*(-1 + \[ScriptM])*\[Omega]*((3*I)*n^2 + 
                n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
                \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                  M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                     \[Omega]^2))))*SphericalHarmonicY[\[ScriptL], 
              \[ScriptM], Pi/2, 0] + ang*\[Omega]*((3*I)*n^2 + 
              n*\[Omega]*(-12 + M*(-12 - I*(2 + a0)*\[Omega])) + 
              \[Omega]*(27 - 6*(2 + a0)*M^2*\[Omega]^2 + 
                M*(27 + (39*I)*(2 + a0)*\[Omega] + 2*(-2 + 5*a0 + 3*a0^2)*
                   \[Omega]^2)))*Conjugate[1/Sqrt[Gamma[-1 + \[ScriptL] - 
                  \[ScriptM]]]]*Sqrt[(2 + \[ScriptL]^2 + 3*\[ScriptM] + 
                \[ScriptM]^2 + \[ScriptL]*(3 + 2*\[ScriptM]))*(\[ScriptL] - 
                 \[ScriptM])!]*SphericalHarmonicY[\[ScriptL], 2 + \[ScriptM], 
              Pi/2, 0]))/(E^((r - rpart)^2/(2*\[Sigma]^2))*enr*n*(1 + n)*
           \[Sigma]*\[Omega]^3)))/(9 - (6*I)*(1 + M)*\[Omega] + 
       (2 + a0)*M*\[Omega]^2)))/((-2 + r)*r^2*
    Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]), 0, 0, 0, 0}
