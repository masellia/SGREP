(* Created with the Wolfram Language : www.wolfram.com *)
{((-2*(2*MBH - r)*(-3*a0^2*r*(-3*MBH + r + n*r) + a0^3*(3*MBH - (1 + n)*r) - 
       r*(r^2*(-3*MBH + r + n*r) + M*(-16*MBH^2 + 14*MBH*r - 3*r^2)) + 
       a0*(-3*r^2*(-3*MBH + r + n*r) + 2*M*(6*MBH^2 - 5*MBH*r + r^2))))/
     (E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
        (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*r^4*
      (a0 + r)^3) + \[Omega]^2)*X[r] + 
  (2*(-2*MBH + r)*(a0^3*MBH + 3*a0^2*MBH*r + MBH*r^3 + 
     M*r*(8*MBH^2 - 6*MBH*r + r^2) + a0*MBH*(4*M*MBH - 2*M*r + 3*r^2))*
    Derivative[1][X][r])/(E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
      (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*r^3*
    (a0 + r)^3) + ((-2*MBH + r)^2*(a0^2 + 4*M*MBH + 2*a0*r - 2*M*r + r^2)*
    Derivative[2][X][r])/(E^(Sqrt[M/(2*a0 - M + 4*MBH)]*
      (Pi - 2*ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r]))*r^2*
    (a0 + r)^2), ((4*I)*ang^2*Sqrt[2*Pi]*(-2*MBH + r)^2*
   (a0^3*(r*(r^2 - r*rpart + 4*\[Sigma]^2) - 
      2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2)) + 
    3*a0^2*r*(r*(r^2 - r*rpart + 4*\[Sigma]^2) - 
      2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2)) + 
    r^3*(r*(r^2 - r*rpart + 4*\[Sigma]^2) - 
      2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2)) - 2*M*(2*MBH - r)*r*
     (-(r*(r^2 - r*rpart + 6*\[Sigma]^2)) + 
      2*MBH*(r^2 - r*rpart + 7*\[Sigma]^2)) + 
    a0*(3*r^2*(r*(r^2 - r*rpart + 4*\[Sigma]^2) - 
        2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2)) - 2*M*(2*MBH - r)*
       (-(r*(r^2 - r*rpart + 5*\[Sigma]^2)) + 
        2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2))))*
   Conjugate[1/Sqrt[(-1 + \[ScriptL])*\[ScriptL]*(1 + \[ScriptL])*
       (2 + \[ScriptL])]]*Conjugate[
    (\[ScriptM]*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]]*
      Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[\[ScriptL], 
       1 + \[ScriptM], Pi/2, 0])/(Sqrt[Gamma[\[ScriptL] - \[ScriptM]]]*
      Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]])])/
  (E^((r^2 - 2*r*rpart + rpart^2 + 4*Sqrt[M/(2*a0 - M + 4*MBH)]*Pi*
       \[Sigma]^2 - 8*Sqrt[M/(2*a0 - M + 4*MBH)]*\[Sigma]^2*
       ArcTan[Sqrt[M*(2*a0 - M + 4*MBH)], a0 - M + r])/(2*\[Sigma]^2))*enr*
   Sqrt[n*(1 + n)]*r^7*(a0 + r)^3*\[Sigma]^3), 
 Xh[0] + (I*(-1 + 2*n)*(-2 + r)*Xh[0])/
   (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
        (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
     \[Omega]) + (-2 + r)^2*
   (((I/8)*ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
         (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[Pi/2]*
      Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/(enr*n*(1 + n)*\[Sigma]*
      (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
       (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])) + ((3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
       (I*(1 + n)*(-1 + 2*n))/
        (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
       (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
        (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*Xh[0])/
     (8*(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) + (-2 + r)^3*
   (\[Mu]*(-1/4*(ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + 
           Sqrt[M/(4 + 2*a0 - M)]*(Pi + 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[Pi/2]*(2 - rpart)*
         \[ScriptM]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
            SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(enr*n*(1 + n)*\[Sigma]^3*
         ((-3*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          4*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
          (6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           4*a0*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           a0^2*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
         I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
          (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
             2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           a0^2*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
             2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                  \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                 Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(48*(2 + a0)^2*enr*n*(1 + n)*
        \[Sigma]*(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
    (((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 12*a0^2*(18 + (-12 + M)*n) - 
          8*a0*(36 + M - 24*n + 7*M*n) + 16*(-9 + 6*n + M^2*n - 
            M*(1 + 4*n))))/(2 + a0)^4 + (12*n*(-1 + 2*n))/
        (-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) - 
       (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + a0^3 + 
          4*M)*(-1 + 2*n)*\[Omega])/((2 + a0)^3*
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/((2 + a0)^3*
         (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
       (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) - 
       ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
         (3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*Xh[0])/
     (64*(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) + (-2 + r)^4*
   (\[Mu]*(((I/32)*ang^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Pi*
        \[ScriptM]*((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
           Sqrt[2*Pi]*\[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
           Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)*
        (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (-2*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
          (84*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + 84*a0*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                 \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] + 
           21*a0^2*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/
              Gamma[\[ScriptL] - \[ScriptM]]] - 
           16*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/
              Gamma[\[ScriptL] - \[ScriptM]]]) + 
         I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4*a0*(-26 + n)*
            Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + a0^2*(-26 + n)*
            Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + 4*(-26*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                   \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] + 
             n*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                  \[ScriptM]]] + 6*Sqrt[(M^3*Gamma[1 + \[ScriptL] - 
                   \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]]))))/
       (64*(2 + a0)^2*enr*Sqrt[M]*n*(1 + n)*\[Sigma]^3*
        (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (((-3*I)*(2 + a0)*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + 
            a0*(72 - 4*(-9 + M)*n) + 8*(6 + 3*M + 3*n - M*n))*
           Sqrt[\[ScriptL] - \[ScriptM]])/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
           (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
         (12*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + 
            a0^3 + 4*M)*Sqrt[\[ScriptL] - \[ScriptM]]*\[Omega])/
          (-E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
           (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
         (6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
           ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
             (6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                    \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                    \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                    \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                  Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*(-9 + 
                Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) - 
            E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
             (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                      \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + a0^2*(9*(-4 + Sqrt[((\[ScriptL] - 
                      \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                     \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
             Gamma[\[ScriptL] - \[ScriptM]]])/
          (-3*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
           (10*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
           8*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
         ((2 + a0)^2*(10 + n)*((-8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (5*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
               2)*\[Omega]*(6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*
                   Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                    \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                    \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                  Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*(-9 + 
                Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
            E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
             (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                      \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + a0^2*(9*(-4 + Sqrt[((\[ScriptL] - 
                      \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                     \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
             Gamma[\[ScriptL] - \[ScriptM]]])/
          ((-3*I)*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
           10*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
           (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
         ((12*I)*(168*a0^3*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                \[ScriptL] - \[ScriptM]]] + 21*a0^4*
             Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                 \[ScriptM]]] + 14*a0*(48*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                    \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 
              7*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                  \[ScriptL] - \[ScriptM]]]) + 24*a0^2*
             (21*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                  \[ScriptL] - \[ScriptM]]] - Sqrt[(M^3*Gamma[
                  1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                  \[ScriptM]]]) + 4*(84*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                    \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 
              25*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/
                 Gamma[\[ScriptL] - \[ScriptM]]] + 4*Sqrt[
                (M^5*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                   \[ScriptM]]])))/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
           Sqrt[M])))/(192*(2 + a0)^4*enr*n*(1 + n)*\[Sigma]*
        (4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
    (((16*(-180*a0^5*(-2 + n) - 15*a0^6*(-2 + n) + 
          36*a0^4*(50 + (-25 + M)*n) + 6*a0^3*(-400*(-2 + n) + 
            9*M*(1 + 6*n)) - 24*a0^2*(150*(-2 + n) + 3*M^2*n - 
            2*M*(8 + 23*n)) + 16*(-60*(-2 + n) + 6*M^3*n + 6*M*(7 + 10*n) - 
            M^2*(5 + 32*n)) - 8*a0*(360*(-2 + n) + M^2*(5 + 50*n) - 
            3*M*(37 + 70*n))))/(2 + a0)^6 + 
       (16*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(120*a0^3 + 30*a0^4 + 
          3*a0^5 + 12*a0^2*(20 + M) + 48*a0*(5 + 2*M) + 16*(6 + 9*M - 2*M^2))*
         (-1 + 2*n)*\[Omega])/((2 + a0)^5*
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((32*I)*(-1 + 2*n)*(30*a0^4*(-2 + 3*n) + a0^5*(-6 + 9*n) - 
          6*a0^3*(40 + 3*(-20 + M)*n) + 12*a0*(-40 + M*(11 - 22*n) + 60*n + 
            2*M^2*n) + 8*(-24 - 24*M*(-1 + n) + 36*n + M^2*(-5 + 6*n)) - 
          6*a0^2*(80 - 120*n + M*(-3 + 20*n))))/((2 + a0)^5*
         (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
       (12*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + a0*(72 - 4*(-9 + M)*n) + 
          8*(6 + 3*M + 3*n - M*n))*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
          (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/((2 + a0)^3*
         (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((48*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + a0^3 + 
          4*M)*\[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
          (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/((2 + a0)^3*
         (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (3*(10 + n)*((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
             12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
             16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
          (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
           ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
           ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) - ((18*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega]*((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
             12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
             16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
          (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
           ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
           ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]))*Xh[0])/
     (768*(4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))), 
 ((-I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*r*(a0 + r)*\[Omega]*Xh[0])/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) + 
  (I*(-1 + 2*n)*Xh[0])/
   (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
        (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
     \[Omega]) + 
  (E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(-1 + 2*n)*r*(a0 + r)*\[Omega]*Xh[0])/
   (Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]*
    (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*
         (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
      \[Omega])) + 2*(-2 + r)*
   (((I/8)*ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
         (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[Pi/2]*
      Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/(enr*n*(1 + n)*\[Sigma]*
      (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
       (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
            2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega])) + ((3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
       (I*(1 + n)*(-1 + 2*n))/
        (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
       (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
        (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*Xh[0])/
     (8*(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) - 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(-2 + r)*r*(a0 + r)*\[Omega]*
    (((I/8)*ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
          (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
       Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/(enr*n*(1 + n)*\[Sigma]*
       (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
        (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
     ((3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
        (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*Xh[0])/
      (8*(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))))/
   Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2] + 
  3*(-2 + r)^2*
   (\[Mu]*(-1/4*(ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + 
           Sqrt[M/(4 + 2*a0 - M)]*(Pi + 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[Pi/2]*(2 - rpart)*
         \[ScriptM]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
            SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(enr*n*(1 + n)*\[Sigma]^3*
         ((-3*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          4*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
          (6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           4*a0*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           a0^2*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                  \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
         I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
          (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
             2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           a0^2*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
             2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
           4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                  \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                 Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(48*(2 + a0)^2*enr*n*(1 + n)*
        \[Sigma]*(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
    (((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 12*a0^2*(18 + (-12 + M)*n) - 
          8*a0*(36 + M - 24*n + 7*M*n) + 16*(-9 + 6*n + M^2*n - 
            M*(1 + 4*n))))/(2 + a0)^4 + (12*n*(-1 + 2*n))/
        (-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) - 
       (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + a0^3 + 
          4*M)*(-1 + 2*n)*\[Omega])/((2 + a0)^3*
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/((2 + a0)^3*
         (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
       (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) - 
       ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
         (3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
           (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
        (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))*Xh[0])/
     (64*(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) - 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(-2 + r)^2*r*(a0 + r)*\[Omega]*
    (\[Mu]*(-1/4*(ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + 
            Sqrt[M/(4 + 2*a0 - M)]*(Pi + 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))*Sqrt[Pi/2]*(2 - rpart)*
          \[ScriptM]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
             SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
            Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
          Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
              \[ScriptM]]])/(enr*n*(1 + n)*\[Sigma]^3*
          ((-3*I)*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
           4*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
            (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         Sqrt[Pi/2]*\[ScriptM]*Conjugate[
          (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
             \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
           (6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
            4*a0*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
            a0^2*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                   \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
          I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
           (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
            a0^2*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                    \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
            4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                   \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                   \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                  Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                   \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
           Gamma[\[ScriptL] - \[ScriptM]]])/(48*(2 + a0)^2*enr*n*(1 + n)*
         \[Sigma]*(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
         (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
     (((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
           12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
           16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
        (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
         ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
         ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
           (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
           (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
           (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]))*Xh[0])/
      (64*(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega]))))/Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2] + 
  4*(-2 + r)^3*
   (\[Mu]*(((I/32)*ang^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
          (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Pi*
        \[ScriptM]*((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
           Sqrt[2*Pi]*\[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
           Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)*
        (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (-2*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
          (84*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + 84*a0*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                 \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] + 
           21*a0^2*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/
              Gamma[\[ScriptL] - \[ScriptM]]] - 
           16*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/
              Gamma[\[ScriptL] - \[ScriptM]]]) + 
         I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4*a0*(-26 + n)*
            Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + a0^2*(-26 + n)*
            Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                \[ScriptM]]] + 4*(-26*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                   \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] + 
             n*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                  \[ScriptM]]] + 6*Sqrt[(M^3*Gamma[1 + \[ScriptL] - 
                   \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]]))))/
       (64*(2 + a0)^2*enr*Sqrt[M]*n*(1 + n)*\[Sigma]^3*
        (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
        (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
         (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
      (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
        Sqrt[Pi/2]*\[ScriptM]*Conjugate[
         (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
            \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        (((-3*I)*(2 + a0)*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + 
            a0*(72 - 4*(-9 + M)*n) + 8*(6 + 3*M + 3*n - M*n))*
           Sqrt[\[ScriptL] - \[ScriptM]])/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
           (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
         (12*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + 
            a0^3 + 4*M)*Sqrt[\[ScriptL] - \[ScriptM]]*\[Omega])/
          (-E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
           (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
         (6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
           ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
             (6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                    \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                    \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                    \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                  Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*(-9 + 
                Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) - 
            E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
             (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                      \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + a0^2*(9*(-4 + Sqrt[((\[ScriptL] - 
                      \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                     \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
             Gamma[\[ScriptL] - \[ScriptM]]])/
          (-3*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
           (10*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
           8*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
         ((2 + a0)^2*(10 + n)*((-8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (5*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
               2)*\[Omega]*(6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*
                   Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                    \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                    \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                  Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*(-9 + 
                Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                     \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
            E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
             (4*a0*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                      \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + a0^2*(9*(-4 + Sqrt[((\[ScriptL] - 
                      \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
              4*(-36 + 6*M + 9*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                     \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]] + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]])))*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
             Gamma[\[ScriptL] - \[ScriptM]]])/
          ((-3*I)*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
           10*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
           (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
         ((12*I)*(168*a0^3*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                \[ScriptL] - \[ScriptM]]] + 21*a0^4*
             Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                 \[ScriptM]]] + 14*a0*(48*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                    \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 
              7*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                  \[ScriptL] - \[ScriptM]]]) + 24*a0^2*
             (21*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                  \[ScriptL] - \[ScriptM]]] - Sqrt[(M^3*Gamma[
                  1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                  \[ScriptM]]]) + 4*(84*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                    \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 
              25*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/
                 Gamma[\[ScriptL] - \[ScriptM]]] + 4*Sqrt[
                (M^5*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                   \[ScriptM]]])))/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
           Sqrt[M])))/(192*(2 + a0)^4*enr*n*(1 + n)*\[Sigma]*
        (4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
    (((16*(-180*a0^5*(-2 + n) - 15*a0^6*(-2 + n) + 
          36*a0^4*(50 + (-25 + M)*n) + 6*a0^3*(-400*(-2 + n) + 
            9*M*(1 + 6*n)) - 24*a0^2*(150*(-2 + n) + 3*M^2*n - 
            2*M*(8 + 23*n)) + 16*(-60*(-2 + n) + 6*M^3*n + 6*M*(7 + 10*n) - 
            M^2*(5 + 32*n)) - 8*a0*(360*(-2 + n) + M^2*(5 + 50*n) - 
            3*M*(37 + 70*n))))/(2 + a0)^6 + 
       (16*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(120*a0^3 + 30*a0^4 + 
          3*a0^5 + 12*a0^2*(20 + M) + 48*a0*(5 + 2*M) + 16*(6 + 9*M - 2*M^2))*
         (-1 + 2*n)*\[Omega])/((2 + a0)^5*
         (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((32*I)*(-1 + 2*n)*(30*a0^4*(-2 + 3*n) + a0^5*(-6 + 9*n) - 
          6*a0^3*(40 + 3*(-20 + M)*n) + 12*a0*(-40 + M*(11 - 22*n) + 60*n + 
            2*M^2*n) + 8*(-24 - 24*M*(-1 + n) + 36*n + M^2*(-5 + 6*n)) - 
          6*a0^2*(80 - 120*n + M*(-3 + 20*n))))/((2 + a0)^5*
         (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
       (12*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + a0*(72 - 4*(-9 + M)*n) + 
          8*(6 + 3*M + 3*n - M*n))*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
          (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/((2 + a0)^3*
         (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       ((48*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + a0^3 + 
          4*M)*\[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
          (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/((2 + a0)^3*
         (1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (3*(10 + n)*((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
             12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
             16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
          (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
           ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
           ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]) - ((18*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega]*((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
             12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
             16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
          (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
           ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
           ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
             (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(
                -1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                   (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                  2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]))*Xh[0])/
     (768*(4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
           (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
        \[Omega]))) - 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(-2 + r)^3*r*(a0 + r)*\[Omega]*
    (\[Mu]*(((I/32)*ang^2*E^(Sqrt[M/(4 + 2*a0 - M)]*
           (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*Pi*
         \[ScriptM]*((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
            Sqrt[2*Pi]*\[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
            Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
          (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
             \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(enr*n*(1 + n)*(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
       (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
            (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*Conjugate[
          (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
             \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         (-2*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]*
           (84*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                \[ScriptL] - \[ScriptM]]] + 84*a0*
             Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                 \[ScriptM]]] + 21*a0^2*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                  \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 
            16*Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                \[ScriptL] - \[ScriptM]]]) + 
          I*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4*a0*(-26 + n)*
             Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                 \[ScriptM]]] + a0^2*(-26 + n)*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                  \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] + 
            4*(-26*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/
                 Gamma[\[ScriptL] - \[ScriptM]]] + n*Sqrt[
                (M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                   \[ScriptM]]] + 6*Sqrt[(M^3*Gamma[1 + \[ScriptL] - 
                    \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]]))))/
        (64*(2 + a0)^2*enr*Sqrt[M]*n*(1 + n)*\[Sigma]^3*
         (E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          I*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])*
         (3*E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
          (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) + 
       (ang^2*E^(-1/2*(2 - rpart)^2/\[Sigma]^2 + Sqrt[M/(4 + 2*a0 - M)]*
            (Pi + 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))*
         Sqrt[Pi/2]*\[ScriptM]*Conjugate[
          (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
             \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
           Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
         (((-3*I)*(2 + a0)*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + 
             a0*(72 - 4*(-9 + M)*n) + 8*(6 + 3*M + 3*n - M*n))*
            Sqrt[\[ScriptL] - \[ScriptM]])/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
            (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (12*(2 + a0)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(8 + 12*a0 + 6*a0^2 + 
             a0^3 + 4*M)*Sqrt[\[ScriptL] - \[ScriptM]]*\[Omega])/
           (-E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
            (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
          (6*(2 + a0)^2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega]*((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 
                  2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]*(6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                     \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                   Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*
                (-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) - 
             E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4*a0*
                (9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                        \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                       \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
               a0^2*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                       \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                       \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                     Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 4*(-36 + 6*M + 9*Sqrt[
                   ((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]] + 
                 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                       \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])))*
            Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]])/(-3*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) + 
            (10*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
            8*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
          ((2 + a0)^2*(10 + n)*((-8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (5*Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                2)*\[Omega]*(6*M + 4*(-9 + Sqrt[((\[ScriptL] - \[ScriptM])*
                    Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                     \[ScriptM]]]) + 4*a0*(-9 + Sqrt[((\[ScriptL] - 
                     \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                   Gamma[1 + \[ScriptL] - \[ScriptM]]]) + a0^2*
                (-9 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                      \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])) + 
             E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*(4*a0*
                (9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                        \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
                 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                       \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]]) + 
               a0^2*(9*(-4 + Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[
                       \[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                       \[ScriptM]]]) + 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*
                     Gamma[\[ScriptL] - \[ScriptM]])/Gamma[1 + \[ScriptL] - 
                      \[ScriptM]]]) + 4*(-36 + 6*M + 9*Sqrt[
                   ((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - \[ScriptM]])/
                    Gamma[1 + \[ScriptL] - \[ScriptM]]] + 
                 2*n*Sqrt[((\[ScriptL] - \[ScriptM])*Gamma[\[ScriptL] - 
                       \[ScriptM]])/Gamma[1 + \[ScriptL] - \[ScriptM]]])))*
            Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]])/((-3*I)*E^(4*Sqrt[M/(4 + 2*a0 - M)]*Pi) - 
            10*E^((Sqrt[M/(4 + 2*a0 - M)]*(9*Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega] + 
            (8*I)*E^(Sqrt[M/(4 + 2*a0 - M)]*(5*Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))*\[Omega]^2) + 
          ((12*I)*(168*a0^3*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/
                Gamma[\[ScriptL] - \[ScriptM]]] + 21*a0^4*
              Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                  \[ScriptM]]] + 14*a0*(48*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                     \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 7*
                Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                   \[ScriptL] - \[ScriptM]]]) + 24*a0^2*
              (21*Sqrt[(M*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                   \[ScriptL] - \[ScriptM]]] - Sqrt[
                (M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                   \[ScriptM]]]) + 4*(84*Sqrt[(M*Gamma[1 + \[ScriptL] - 
                     \[ScriptM]])/Gamma[\[ScriptL] - \[ScriptM]]] - 25*
                Sqrt[(M^3*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[
                   \[ScriptL] - \[ScriptM]]] + 4*Sqrt[
                 (M^5*Gamma[1 + \[ScriptL] - \[ScriptM]])/Gamma[\[ScriptL] - 
                    \[ScriptM]]])))/(E^(2*Sqrt[M/(4 + 2*a0 - M)]*Pi)*
            Sqrt[M])))/(192*(2 + a0)^4*enr*n*(1 + n)*\[Sigma]*
         (4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                  Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]))) + 
     (((16*(-180*a0^5*(-2 + n) - 15*a0^6*(-2 + n) + 
           36*a0^4*(50 + (-25 + M)*n) + 6*a0^3*(-400*(-2 + n) + 
             9*M*(1 + 6*n)) - 24*a0^2*(150*(-2 + n) + 3*M^2*n - 
             2*M*(8 + 23*n)) + 16*(-60*(-2 + n) + 6*M^3*n + 6*M*(7 + 10*n) - 
             M^2*(5 + 32*n)) - 8*a0*(360*(-2 + n) + M^2*(5 + 50*n) - 
             3*M*(37 + 70*n))))/(2 + a0)^6 + 
        (16*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                 Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(120*a0^3 + 30*a0^4 + 
           3*a0^5 + 12*a0^2*(20 + M) + 48*a0*(5 + 2*M) + 
           16*(6 + 9*M - 2*M^2))*(-1 + 2*n)*\[Omega])/
         ((2 + a0)^5*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega])) + ((32*I)*(-1 + 2*n)*(30*a0^4*(-2 + 3*n) + 
           a0^5*(-6 + 9*n) - 6*a0^3*(40 + 3*(-20 + M)*n) + 
           12*a0*(-40 + M*(11 - 22*n) + 60*n + 2*M^2*n) + 
           8*(-24 - 24*M*(-1 + n) + 36*n + M^2*(-5 + 6*n)) - 
           6*a0^2*(80 - 120*n + M*(-3 + 20*n))))/((2 + a0)^5*
          (2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                   Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])) - 
        (12*(18*a0^2*(2 + n) + 3*a0^3*(2 + n) + a0*(72 - 4*(-9 + M)*n) + 
           8*(6 + 3*M + 3*n - M*n))*(3 - ((3*(2 + a0)^2 - 4*M)*n)/
            (2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
            (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
           (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
            (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
         ((2 + a0)^3*(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega])) + ((48*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*\[Omega]*
          (3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + (I*(1 + n)*(-1 + 2*n))/
            (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega]) + 
           (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                    Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*(-1 + 2*n)*\[Omega])/
            (I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[(2 + a0 - M)/
                     Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*\[Omega])))/
         ((2 + a0)^3*(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
            \[Omega])) + (3*(10 + n)*
          ((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 12*a0^2*(18 + 
                (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
              16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
           (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
            ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
            ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]) - ((18*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
             (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
          \[Omega]*((4*(24*a0^3*(-3 + 2*n) + a0^4*(-9 + 6*n) - 
              12*a0^2*(18 + (-12 + M)*n) - 8*a0*(36 + M - 24*n + 7*M*n) + 
              16*(-9 + 6*n + M^2*n - M*(1 + 4*n))))/(2 + a0)^4 + 
           (12*n*(-1 + 2*n))/(-1 + (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) - (8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             (8 + 12*a0 + 6*a0^2 + a0^3 + 4*M)*(-1 + 2*n)*\[Omega])/
            ((2 + a0)^3*(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega])) + ((32*I)*M*(-1 + 2*n)*(-2 + (2 + a0)*n))/
            ((2 + a0)^3*(2*I + 8*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                   2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
               \[Omega])) + (9*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) + (2*n*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega]) - ((8*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
             \[Omega]*(3 - ((3*(2 + a0)^2 - 4*M)*n)/(2 + a0)^2 + 
              (I*(1 + n)*(-1 + 2*n))/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega]) + (2*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 
                    2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
                (-1 + 2*n)*\[Omega])/(I + 4*E^((Sqrt[M/(4 + 2*a0 - M)]*
                    (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/
                   2)*\[Omega])))/(1 - (2*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
                 (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
              \[Omega])))/(9/4 - (3*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
              (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
           \[Omega]))*Xh[0])/
      (768*(4 - (4*I)*E^((Sqrt[M/(4 + 2*a0 - M)]*
            (Pi - 2*ArcTan[(2 + a0 - M)/Sqrt[-(M*(-4 - 2*a0 + M))]]))/2)*
         \[Omega]))))/Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2], 
 X\[Infinity][0] + (I*(1 + n)*X\[Infinity][0])/(r*\[Omega]) + 
  ((ang^2*Sqrt[Pi/2]*(r - rpart)*\[ScriptM]*\[Mu]*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
      Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
        Gamma[\[ScriptL] - \[ScriptM]]])/(E^((r - rpart)^2/(2*\[Sigma]^2))*
      \[Sigma]^3*(enr*n*\[Omega] + enr*n^2*\[Omega])) - 
    ((n + n^2 + (3*I)*(1 + M)*\[Omega])*X\[Infinity][0])/(2*\[Omega]^2))/
   r^2 + ((ang^2*Pi*\[ScriptM]*\[Mu]*((4*Sqrt[2/Pi]*\[Omega])/
        (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
       ((r - rpart)*(2*I - I*n + 8*(1 + M)*\[Omega]))/
        (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
      Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
        Gamma[\[ScriptL] - \[ScriptM]]])/(3*enr*n*(1 + n)*\[Omega]^2) - 
    ((I/6)*(-n^2 + n^3 + n*(-2 + (3*I)*(1 + M)*\[Omega] + 
         2*(2 + a0)*M*\[Omega]^2) - 6*\[Omega]*
        (I + M*(I + 2*(2 + a0)*\[Omega])))*X\[Infinity][0])/\[Omega]^3)/r^3 + 
  (-1/12*(ang^2*Pi*\[ScriptM]*\[Mu]*
       ((4*Sqrt[2/Pi]*\[Omega]*(5*I - I*n + 15*(1 + M)*\[Omega]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(-10 - n^2 + (55*I)*(1 + M)*\[Omega] + 
           6*(4 + 18*M + 5*a0*M + 4*M^2)*\[Omega]^2 + 
           n*(7 - (8*I)*(1 + M)*\[Omega])))/(E^((r - rpart)^2/(2*\[Sigma]^2))*
          Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
        (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
           \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
       Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
         Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)*\[Omega]^3) + 
    ((-6*n^3 + n^4 + n^2*(3 - (12*I)*(1 + M)*\[Omega] + 
         8*(2 + a0)*M*\[Omega]^2) + 2*n*(5 - (18*I)*(1 + M)*\[Omega] - 
         8*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
          \[Omega]^3) - (5*I)*\[Omega]*(-6 + (9*I)*\[Omega] + 
         M^2*\[Omega]*(9*I - 2*(2 + a0)*\[Omega]) + 
         6*M*(-1 + I*(7 + 2*a0)*\[Omega] + (2 + 5*a0 + 2*a0^2)*\[Omega]^2)))*
      X\[Infinity][0])/(24*\[Omega]^4))/r^4, 
 -((E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
           a0 - M + r]))/2)*(1 + n)*(a0 + r)*X\[Infinity][0])/
    ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2])) - 
  (I*(1 + n)*X\[Infinity][0])/(r^2*\[Omega]) + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*r*(a0 + r)*\[Omega]*X\[Infinity][0])/
   ((-2 + r)*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (2*((ang^2*Sqrt[Pi/2]*(r - rpart)*\[ScriptM]*\[Mu]*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
       Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
         Gamma[\[ScriptL] - \[ScriptM]]])/(E^((r - rpart)^2/(2*\[Sigma]^2))*
       \[Sigma]^3*(enr*n*\[Omega] + enr*n^2*\[Omega])) - 
     ((n + n^2 + (3*I)*(1 + M)*\[Omega])*X\[Infinity][0])/(2*\[Omega]^2)))/
   r^3 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(a0 + r)*\[Omega]*
    ((ang^2*Sqrt[Pi/2]*(r - rpart)*\[ScriptM]*\[Mu]*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
       Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
         Gamma[\[ScriptL] - \[ScriptM]]])/(E^((r - rpart)^2/(2*\[Sigma]^2))*
       \[Sigma]^3*(enr*n*\[Omega] + enr*n^2*\[Omega])) - 
     ((n + n^2 + (3*I)*(1 + M)*\[Omega])*X\[Infinity][0])/(2*\[Omega]^2)))/
   ((-2 + r)*r*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (3*((ang^2*Pi*\[ScriptM]*\[Mu]*((4*Sqrt[2/Pi]*\[Omega])/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(2*I - I*n + 8*(1 + M)*\[Omega]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
       Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
         Gamma[\[ScriptL] - \[ScriptM]]])/(3*enr*n*(1 + n)*\[Omega]^2) - 
     ((I/6)*(-n^2 + n^3 + n*(-2 + (3*I)*(1 + M)*\[Omega] + 
          2*(2 + a0)*M*\[Omega]^2) - 6*\[Omega]*
         (I + M*(I + 2*(2 + a0)*\[Omega])))*X\[Infinity][0])/\[Omega]^3))/
   r^4 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(a0 + r)*\[Omega]*
    ((ang^2*Pi*\[ScriptM]*\[Mu]*((4*Sqrt[2/Pi]*\[Omega])/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(2*I - I*n + 8*(1 + M)*\[Omega]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
       Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
         Gamma[\[ScriptL] - \[ScriptM]]])/(3*enr*n*(1 + n)*\[Omega]^2) - 
     ((I/6)*(-n^2 + n^3 + n*(-2 + (3*I)*(1 + M)*\[Omega] + 
          2*(2 + a0)*M*\[Omega]^2) - 6*\[Omega]*
         (I + M*(I + 2*(2 + a0)*\[Omega])))*X\[Infinity][0])/\[Omega]^3))/
   ((-2 + r)*r^2*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2]) - 
  (4*(-1/12*(ang^2*Pi*\[ScriptM]*\[Mu]*
        ((4*Sqrt[2/Pi]*\[Omega]*(5*I - I*n + 15*(1 + M)*\[Omega]))/
          (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
         ((r - rpart)*(-10 - n^2 + (55*I)*(1 + M)*\[Omega] + 
            6*(4 + 18*M + 5*a0*M + 4*M^2)*\[Omega]^2 + 
            n*(7 - (8*I)*(1 + M)*\[Omega])))/
          (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)*\[Omega]^3) + 
     ((-6*n^3 + n^4 + n^2*(3 - (12*I)*(1 + M)*\[Omega] + 
          8*(2 + a0)*M*\[Omega]^2) + 2*n*(5 - (18*I)*(1 + M)*\[Omega] - 
          8*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
           \[Omega]^3) - (5*I)*\[Omega]*(-6 + (9*I)*\[Omega] + 
          M^2*\[Omega]*(9*I - 2*(2 + a0)*\[Omega]) + 
          6*M*(-1 + I*(7 + 2*a0)*\[Omega] + (2 + 5*a0 + 2*a0^2)*\[Omega]^2)))*
       X\[Infinity][0])/(24*\[Omega]^4)))/r^5 + 
  (I*E^((Sqrt[M/(4 + 2*a0 - M)]*(Pi - 2*ArcTan[Sqrt[(4 + 2*a0 - M)*M], 
          a0 - M + r]))/2)*(a0 + r)*\[Omega]*
    (-1/12*(ang^2*Pi*\[ScriptM]*\[Mu]*
        ((4*Sqrt[2/Pi]*\[Omega]*(5*I - I*n + 15*(1 + M)*\[Omega]))/
          (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
         ((r - rpart)*(-10 - n^2 + (55*I)*(1 + M)*\[Omega] + 
            6*(4 + 18*M + 5*a0*M + 4*M^2)*\[Omega]^2 + 
            n*(7 - (8*I)*(1 + M)*\[Omega])))/
          (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)*\[Omega]^3) + 
     ((-6*n^3 + n^4 + n^2*(3 - (12*I)*(1 + M)*\[Omega] + 
          8*(2 + a0)*M*\[Omega]^2) + 2*n*(5 - (18*I)*(1 + M)*\[Omega] - 
          8*(2 + a0)*M*\[Omega]^2 + (2*I)*(2 + a0)*(3*a0 - 2*M)*M*
           \[Omega]^3) - (5*I)*\[Omega]*(-6 + (9*I)*\[Omega] + 
          M^2*\[Omega]*(9*I - 2*(2 + a0)*\[Omega]) + 
          6*M*(-1 + I*(7 + 2*a0)*\[Omega] + (2 + 5*a0 + 2*a0^2)*\[Omega]^2)))*
       X\[Infinity][0])/(24*\[Omega]^4)))/
   ((-2 + r)*r^3*Sqrt[a0^2 - 2*M*(-2 + r) + 2*a0*r + r^2])}
