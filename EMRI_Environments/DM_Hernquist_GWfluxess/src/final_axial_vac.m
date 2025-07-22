(* Created with the Wolfram Language : www.wolfram.com *)
{((2*(2*MBH - r)*(-3*MBH + r + n*r))/r^4 + \[Omega]^2)*X[r] + 
  (2*MBH*(-2*MBH + r)*Derivative[1][X][r])/r^3 + 
  ((-2*MBH + r)^2*Derivative[2][X][r])/r^2, 
 ((4*I)*ang^2*Sqrt[2*Pi]*(-2*MBH + r)^2*(r*(r^2 - r*rpart + 4*\[Sigma]^2) - 
    2*MBH*(r^2 - r*rpart + 6*\[Sigma]^2))*
   Conjugate[1/Sqrt[(-1 + \[ScriptL])*\[ScriptL]*(1 + \[ScriptL])*
       (2 + \[ScriptL])]]*Conjugate[
    (\[ScriptM]*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]]*
      Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[\[ScriptL], 
       1 + \[ScriptM], Pi/2, 0])/(Sqrt[Gamma[\[ScriptL] - \[ScriptM]]]*
      Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]])])/
  (E^((r^2 - 2*r*rpart + rpart^2)/(2*\[Sigma]^2))*enr*Sqrt[n*(1 + n)]*r^7*
   \[Sigma]^3), Xh[0] + (I*(-1 + 2*n)*(-2 + r)*Xh[0])/(2*I + 8*\[Omega]) + 
  (-2 + r)^2*(-1/8*(ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
       \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(n + n^2)*\[Sigma]*
       (I + 2*\[Omega])) - ((1 + n*(-1 + n + (4*I)*\[Omega]) - 
       (5*I)*\[Omega])*Xh[0])/(4*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2))) + 
  (-2 + r)^3*((ang^2*Pi*\[ScriptM]*\[Mu]*Conjugate[
       (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
          \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
      ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
       ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
          (Sqrt[\[ScriptL] - \[ScriptM]] - 
           9*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]) - (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
           4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
         \[Sigma])))/(48*enr*n*(1 + n)*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
    ((I/24)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 
       96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2))*Xh[0])/
     ((I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega]))) + 
  (-2 + r)^4*(((3*ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
        \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (32*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (9*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (3*ang^2*n^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + ((I/8)*ang^2*Pi*\[ScriptM]*\[Mu]*
        ((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^3))*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)) + 
      (((3*I)/4)*ang^2*Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]^3) + (((21*I)/16)*ang^2*Sqrt[Pi/2]*\[ScriptM]*
        \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]) + (11*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (5*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(96*enr*n*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (ang^2*n*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) - ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
      ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*n*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])))/
     (4/(1 + n) + (4*n)/(1 + n) - ((4*I)*\[Omega])/(1 + n) - 
      ((4*I)*n*\[Omega])/(1 + n)) + 
    ((5/(8*(1 + n)) + (5*n)/(16*(1 + n)) - (5*n^2)/(16*(1 + n)) - 
       ((I/4)*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((I/8)*n*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       (((3*I)/8)*n^2*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (n*(-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (3*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (16*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (9*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (3*n^2*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       (((5*I)/48)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (((11*I)/96)*n*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        ((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       ((I/96)*n^2*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (\[Omega]*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        (16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       (n*\[Omega]*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/(16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])))*Xh[0])/(4/(1 + n) + (4*n)/(1 + n) - 
      ((4*I)*\[Omega])/(1 + n) - ((4*I)*n*\[Omega])/(1 + n))), 
 ((-I)*r*\[Omega]*Xh[0])/(-2 + r) + (I*(-1 + 2*n)*Xh[0])/(2*I + 8*\[Omega]) + 
  ((-1 + 2*n)*r*\[Omega]*Xh[0])/(2*I + 8*\[Omega]) + 
  2*(-2 + r)*(-1/8*(ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
       \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(n + n^2)*\[Sigma]*
       (I + 2*\[Omega])) - ((1 + n*(-1 + n + (4*I)*\[Omega]) - 
       (5*I)*\[Omega])*Xh[0])/(4*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2))) - 
  I*(-2 + r)*r*\[Omega]*(-1/8*(ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*
       \[ScriptM]*\[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(n + n^2)*\[Sigma]*
       (I + 2*\[Omega])) - ((1 + n*(-1 + n + (4*I)*\[Omega]) - 
       (5*I)*\[Omega])*Xh[0])/(4*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2))) + 
  3*(-2 + r)^2*((ang^2*Pi*\[ScriptM]*\[Mu]*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
      ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
       ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
          (Sqrt[\[ScriptL] - \[ScriptM]] - 
           9*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]) - (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
           4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
         \[Sigma])))/(48*enr*n*(1 + n)*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
    ((I/24)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 
       96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2))*Xh[0])/
     ((I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega]))) - 
  I*(-2 + r)^2*r*\[Omega]*
   ((ang^2*Pi*\[ScriptM]*\[Mu]*Conjugate[
       (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
          \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
      ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
         Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
             \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
       ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
          (Sqrt[\[ScriptL] - \[ScriptM]] - 
           9*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]) - (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
           4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
         \[Sigma])))/(48*enr*n*(1 + n)*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
    ((I/24)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 
       96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2))*Xh[0])/
     ((I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega]))) + 
  4*(-2 + r)^3*(((3*ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
        \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (32*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (9*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (3*ang^2*n^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + ((I/8)*ang^2*Pi*\[ScriptM]*\[Mu]*
        ((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^3))*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)) + 
      (((3*I)/4)*ang^2*Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]^3) + (((21*I)/16)*ang^2*Sqrt[Pi/2]*\[ScriptM]*
        \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]) + (11*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (5*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(96*enr*n*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (ang^2*n*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) - ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
      ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*n*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])))/
     (4/(1 + n) + (4*n)/(1 + n) - ((4*I)*\[Omega])/(1 + n) - 
      ((4*I)*n*\[Omega])/(1 + n)) + 
    ((5/(8*(1 + n)) + (5*n)/(16*(1 + n)) - (5*n^2)/(16*(1 + n)) - 
       ((I/4)*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((I/8)*n*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       (((3*I)/8)*n^2*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (n*(-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (3*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (16*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (9*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (3*n^2*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       (((5*I)/48)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (((11*I)/96)*n*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        ((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       ((I/96)*n^2*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (\[Omega]*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        (16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       (n*\[Omega]*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/(16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])))*Xh[0])/(4/(1 + n) + (4*n)/(1 + n) - 
      ((4*I)*\[Omega])/(1 + n) - ((4*I)*n*\[Omega])/(1 + n))) - 
  I*(-2 + r)^3*r*\[Omega]*
   (((3*ang^2*Sqrt[Pi/2]*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (32*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (9*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + (3*ang^2*n^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (64*E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) - ((I/16)*ang^2*n*Sqrt[Pi/2]*
        Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
       (E^((2 - rpart)^2/(2*\[Sigma]^2))*enr*(1 + n)*(n + n^2)*\[Sigma]*
        (I + 2*\[Omega])) + ((I/8)*ang^2*Pi*\[ScriptM]*\[Mu]*
        ((2 - rpart)^2/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^5) - 1/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma]^3))*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(enr*n*(1 + n)) + 
      (((3*I)/4)*ang^2*Sqrt[Pi/2]*(2 - rpart)*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]^3) + (((21*I)/16)*ang^2*Sqrt[Pi/2]*\[ScriptM]*
        \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/
          Gamma[\[ScriptL] - \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*
        enr*n*(1 + n)*\[Sigma]) + (11*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (5*ang^2*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(96*enr*n*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) + (ang^2*n*Pi*\[ScriptM]*\[Mu]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(192*enr*(1 + n)^2*(I + 2*\[Omega])*
        (3*I + 4*\[Omega])) - ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])) - 
      ((I/32)*ang^2*Pi*\[ScriptM]*\[Mu]*\[Omega]*
        Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
           SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]]*
        ((6*Sqrt[2/Pi]*(2 - rpart)*(I + 2*\[Omega])*
           Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
               \[ScriptM]]])/(E^((2 - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3) + 
         ((-2*I)*n*Sqrt[\[ScriptL] - \[ScriptM]] - 8*\[Omega]*
            (Sqrt[\[ScriptL] - \[ScriptM]] - 9*Sqrt[Gamma[1 + \[ScriptL] - 
                  \[ScriptM]]/Gamma[\[ScriptL] - \[ScriptM]]]) - 
           (9*I)*(Sqrt[\[ScriptL] - \[ScriptM]] - 
             4*Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]/Gamma[\[ScriptL] - 
                  \[ScriptM]]]))/(E^((2 - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*
           \[Sigma])))/(enr*n*(1 + n)^2*(I + 2*\[Omega])*(3*I + 4*\[Omega])))/
     (4/(1 + n) + (4*n)/(1 + n) - ((4*I)*\[Omega])/(1 + n) - 
      ((4*I)*n*\[Omega])/(1 + n)) + 
    ((5/(8*(1 + n)) + (5*n)/(16*(1 + n)) - (5*n^2)/(16*(1 + n)) - 
       ((I/4)*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((I/8)*n*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       (((3*I)/8)*n^2*(-1 + 2*n))/((1 + n)*(2*I + 8*\[Omega])) + 
       ((-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (n*(-1 + 2*n)*\[Omega])/(8*(1 + n)*(2*I + 8*\[Omega])) + 
       (3*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (16*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (9*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) + 
       (3*n^2*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega]))/
        (32*(1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       ((I/8)*n*(1 + n*(-1 + n + (4*I)*\[Omega]) - (5*I)*\[Omega])*\[Omega])/
        ((1 + n)*(-1 + (6*I)*\[Omega] + 8*\[Omega]^2)) - 
       (((5*I)/48)*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (((11*I)/96)*n*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        ((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       ((I/96)*n^2*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/((1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])) - (\[Omega]*(-9 + 2*n^3 + 
          n^2*(-5 + (24*I)*\[Omega]) + (51*I)*\[Omega] + 96*\[Omega]^2 + 
          n*(11 - (42*I)*\[Omega] - 48*\[Omega]^2)))/
        (16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*(3*I + 4*\[Omega])) - 
       (n*\[Omega]*(-9 + 2*n^3 + n^2*(-5 + (24*I)*\[Omega]) + 
          (51*I)*\[Omega] + 96*\[Omega]^2 + n*(11 - (42*I)*\[Omega] - 
            48*\[Omega]^2)))/(16*(1 + n)*(I + 2*\[Omega])*(I + 4*\[Omega])*
         (3*I + 4*\[Omega])))*Xh[0])/(4/(1 + n) + (4*n)/(1 + n) - 
      ((4*I)*\[Omega])/(1 + n) - ((4*I)*n*\[Omega])/(1 + n))), 
 X\[Infinity][0] + (I*(1 + n)*X\[Infinity][0])/(r*\[Omega]) + 
  ((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
      (((4*I)*Sqrt[2/Pi]*(-5 + n + (15*I)*\[Omega])*\[Omega])/
        (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
       ((r - rpart)*(10 + n^2 + n*(-7 + (8*I)*\[Omega]) - 
          \[Omega]*(55*I + 24*\[Omega])))/(E^((r - rpart)^2/(2*\[Sigma]^2))*
         Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
       (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
          \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
     (12*enr*n*(1 + n)*\[Omega]^3) + ((n + n^2 + (3*I)*\[Omega])*
      (10 - 7*n + n^2 - (15*I)*\[Omega])*X\[Infinity][0])/(24*\[Omega]^4))/
   r^4 + ((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
      ((4*Sqrt[2/Pi]*\[Omega])/(E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
       ((r - rpart)*(2*I - I*n + 8*\[Omega]))/
        (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
      Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
     (3*enr*n*(1 + n)*\[Omega]^2) - 
    ((I/6)*(-2 + n)*(n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/\[Omega]^3)/
   r^3 + ((ang^2*Sqrt[Pi/2]*(r - rpart)*Sqrt[\[ScriptL] - \[ScriptM]]*
      \[ScriptM]*\[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
         SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
        Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
     (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3*(enr*n*\[Omega] + 
       enr*n^2*\[Omega])) - ((n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/
     (2*\[Omega]^2))/r^2, -(((1 + n)*X\[Infinity][0])/(-2 + r)) - 
  (I*(1 + n)*X\[Infinity][0])/(r^2*\[Omega]) + (I*r*\[Omega]*X\[Infinity][0])/
   (-2 + r) - (4*((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
       (((4*I)*Sqrt[2/Pi]*(-5 + n + (15*I)*\[Omega])*\[Omega])/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(10 + n^2 + n*(-7 + (8*I)*\[Omega]) - 
           \[Omega]*(55*I + 24*\[Omega])))/(E^((r - rpart)^2/(2*\[Sigma]^2))*
          Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
        (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
           \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (12*enr*n*(1 + n)*\[Omega]^3) + ((n + n^2 + (3*I)*\[Omega])*
       (10 - 7*n + n^2 - (15*I)*\[Omega])*X\[Infinity][0])/(24*\[Omega]^4)))/
   r^5 + (I*\[Omega]*((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
       \[Mu]*(((4*I)*Sqrt[2/Pi]*(-5 + n + (15*I)*\[Omega])*\[Omega])/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(10 + n^2 + n*(-7 + (8*I)*\[Omega]) - 
           \[Omega]*(55*I + 24*\[Omega])))/(E^((r - rpart)^2/(2*\[Sigma]^2))*
          Sqrt[2*Pi]*\[Sigma]^3))*Conjugate[
        (Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
           \[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (12*enr*n*(1 + n)*\[Omega]^3) + ((n + n^2 + (3*I)*\[Omega])*
       (10 - 7*n + n^2 - (15*I)*\[Omega])*X\[Infinity][0])/(24*\[Omega]^4)))/
   ((-2 + r)*r^3) - 
  (3*((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*\[Mu]*
       ((4*Sqrt[2/Pi]*\[Omega])/(E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]) - 
        ((r - rpart)*(2*I - I*n + 8*\[Omega]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (3*enr*n*(1 + n)*\[Omega]^2) - 
     ((I/6)*(-2 + n)*(n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/\[Omega]^3))/
   r^4 + (I*\[Omega]*((ang^2*Pi*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
       \[Mu]*((4*Sqrt[2/Pi]*\[Omega])/(E^((r - rpart)^2/(2*\[Sigma]^2))*
          \[Sigma]) - ((r - rpart)*(2*I - I*n + 8*\[Omega]))/
         (E^((r - rpart)^2/(2*\[Sigma]^2))*Sqrt[2*Pi]*\[Sigma]^3))*
       Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (3*enr*n*(1 + n)*\[Omega]^2) - 
     ((I/6)*(-2 + n)*(n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/\[Omega]^3))/
   ((-2 + r)*r^2) - 
  (2*((ang^2*Sqrt[Pi/2]*(r - rpart)*Sqrt[\[ScriptL] - \[ScriptM]]*\[ScriptM]*
       \[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3*(enr*n*\[Omega] + 
        enr*n^2*\[Omega])) - ((n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/
      (2*\[Omega]^2)))/r^3 + 
  (I*\[Omega]*((ang^2*Sqrt[Pi/2]*(r - rpart)*Sqrt[\[ScriptL] - \[ScriptM]]*
       \[ScriptM]*\[Mu]*Conjugate[(Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*
          SphericalHarmonicY[\[ScriptL], 1 + \[ScriptM], Pi/2, 0])/
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]])/
      (E^((r - rpart)^2/(2*\[Sigma]^2))*\[Sigma]^3*(enr*n*\[Omega] + 
        enr*n^2*\[Omega])) - ((n + n^2 + (3*I)*\[Omega])*X\[Infinity][0])/
      (2*\[Omega]^2)))/((-2 + r)*r)}
