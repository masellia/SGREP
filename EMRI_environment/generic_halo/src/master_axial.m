(* Created with the Wolfram Language : www.wolfram.com *)
{\[Psi][r]*(\[Omega]^2 - 
    (F[r]*(-6*m[r] + r*(\[ScriptL] + \[ScriptL]^2 + Derivative[1][m][r])))/
     r^3) + (F[r]*(2*m[r] - r*Derivative[1][m][r])*Derivative[1][\[Psi]][r])/
   r^2 + (F[r]*(r - 2*m[r])*Derivative[2][\[Psi]][r])/r, 
 ((-8*I)*Lpart^2*Pi*\[Mu]*(Conjugate[\[ScriptM]*Cot[\[Theta]]*
       SphericalHarmonicY[\[ScriptL], \[ScriptM], \[Theta], \[Phi]]] - 
     Conjugate[\[ScriptM]*(\[ScriptM]*Cot[\[Theta]]*SphericalHarmonicY[
          \[ScriptL], \[ScriptM], \[Theta], \[Phi]] + 
        (Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]]*
          Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
           \[ScriptL], 1 + \[ScriptM], \[Theta], \[Phi]])/
         (E^(I*\[Phi])*Sqrt[Gamma[\[ScriptL] - \[ScriptM]]]*
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]))])*Sin[\[Theta]]*
    \[Delta][r]*(4*r*Sqrt[(\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
          \[ScriptL]^3))/r]*F[r]^(3/2)*Sqrt[r - 2*m[r]]*m[r] + 
     (r*Sqrt[r*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3)]*
       (F[r]*(r - 2*m[r]))^(3/2)*Derivative[1][F][r])/F[r] - 
     2*Sqrt[r*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3)]*
      F[r]*Sqrt[F[r]*(r - 2*m[r])]*(-9*m[r] + r*(4 + Derivative[1][m][r]))))/
   (Epart*r*Sqrt[r*(-1 + \[ScriptL])*\[ScriptL]*(1 + \[ScriptL])*
      (2 + \[ScriptL])]*Sqrt[(\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
        \[ScriptL]^3))/r]*
    Sqrt[(r^9*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3))/
      (F[r]*(r - 2*m[r]))]*(r - 2*m[r])) - 
  ((16*I)*Lpart^2*Pi*\[Mu]*(Conjugate[\[ScriptM]*Cot[\[Theta]]*
       SphericalHarmonicY[\[ScriptL], \[ScriptM], \[Theta], \[Phi]]] - 
     Conjugate[\[ScriptM]*(\[ScriptM]*Cot[\[Theta]]*SphericalHarmonicY[
          \[ScriptL], \[ScriptM], \[Theta], \[Phi]] + 
        (Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]]*
          Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
           \[ScriptL], 1 + \[ScriptM], \[Theta], \[Phi]])/
         (E^(I*\[Phi])*Sqrt[Gamma[\[ScriptL] - \[ScriptM]]]*
          Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]))])*F[r]*
    Sqrt[F[r]*(r - 2*m[r])]*Sin[\[Theta]]*Derivative[1][\[Delta]][r])/
   (Epart*Sqrt[(\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3))/
      r]*Sqrt[(r^9*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
        \[ScriptL]^3))/(F[r]*(r - 2*m[r]))]), 
 \[Psi][r]*(\[Omega]^2 - (F[r]*(E^r*\[ScriptL]*(1 + \[ScriptL]) - 6*F[r] + 
       Derivative[1][m][r]))/E^(3*r)) - 
  (F[r]*(E^r - 4*F[r] + Derivative[1][m][r])*Derivative[1][\[Psi]][r])/
   E^(3*r) + ((E^r - 2*F[r])*F[r]*Derivative[2][\[Psi]][r])/E^(3*r), 
 ((-8*I)*Lpart^2*Pi*\[Mu]*(Conjugate[\[ScriptM]*Cot[\[Theta]]*
      SphericalHarmonicY[\[ScriptL], \[ScriptM], \[Theta], \[Phi]]] - 
    Conjugate[\[ScriptM]*(\[ScriptM]*Cot[\[Theta]]*SphericalHarmonicY[
         \[ScriptL], \[ScriptM], \[Theta], \[Phi]] + 
       (Sqrt[Gamma[1 + \[ScriptL] - \[ScriptM]]]*
         Sqrt[Gamma[2 + \[ScriptL] + \[ScriptM]]]*SphericalHarmonicY[
          \[ScriptL], 1 + \[ScriptM], \[Theta], \[Phi]])/
        (E^(I*\[Phi])*Sqrt[Gamma[\[ScriptL] - \[ScriptM]]]*
         Sqrt[Gamma[1 + \[ScriptL] + \[ScriptM]]]))])*Sin[\[Theta]]*
   (4*E^r*Sqrt[(\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3))/
       E^r]*Sqrt[E^r - 2*F[r]]*F[r]^(5/2)*\[Delta][E^r] + 
    E^r*Sqrt[E^r*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
        \[ScriptL]^3)]*Sqrt[(E^r - 2*F[r])*F[r]]*\[Delta][E^r]*
     Derivative[1][F][r] - 
    2*Sqrt[E^r*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3)]*
     F[r]^2*Sqrt[(E^r - 2*F[r])*F[r]]*(-9*\[Delta][E^r] + 
      2*E^r*Derivative[1][\[Delta]][E^r]) + 
    2*Sqrt[E^r*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + \[ScriptL]^3)]*
     F[r]*Sqrt[(E^r - 2*F[r])*F[r]]*
     (-(\[Delta][E^r]*(4*E^r + Derivative[1][F][r] + Derivative[1][m][r])) + 
      E^(2*r)*Derivative[1][\[Delta]][E^r])))/
  (E^r*Epart*Sqrt[E^r*(-1 + \[ScriptL])*\[ScriptL]*(1 + \[ScriptL])*
     (2 + \[ScriptL])]*Sqrt[(\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
       \[ScriptL]^3))/E^r]*(E^r - 2*F[r])*
   Sqrt[(E^(9*r)*\[ScriptL]*(-2 - \[ScriptL] + 2*\[ScriptL]^2 + 
       \[ScriptL]^3))/((E^r - 2*F[r])*F[r])]), 
 -((F[r]*\[Psi][r]*(-6*m[r] + r*(\[ScriptL] + \[ScriptL]^2 + 
        Derivative[1][m][r])))/r^3) + 
  ((2*I)*\[Omega]*Sqrt[(F[r]*(r - 2*m[r]))/r] + 
    (F[r]*(2*m[r] - r*Derivative[1][m][r]))/r^2)*Derivative[1][\[Psi]][r] + 
  (F[r]*(r - 2*m[r])*Derivative[2][\[Psi]][r])/r, 
 -((F[r]*\[Psi][r]*(E^r*\[ScriptL]*(1 + \[ScriptL]) - 6*F[r] + 
      Derivative[1][m][r]))/E^(3*r)) + 
  ((4*F[r]^2 + (2*I)*E^(2*r)*\[Omega]*Sqrt[F[r] - (2*F[r]^2)/E^r] - 
     F[r]*(E^r + Derivative[1][m][r]))*Derivative[1][\[Psi]][r])/E^(3*r) + 
  ((E^r - 2*F[r])*F[r]*Derivative[2][\[Psi]][r])/E^(3*r)}
