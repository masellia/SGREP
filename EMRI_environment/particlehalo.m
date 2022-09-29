(* ::Package:: *)

BeginPackage["particlehalo`"];


(* ::Subsection:: *)
(*Function usage and load*)


waveaxial::usage = 
		"Integrator for axial perturbations - full inhormogeneous integration waveaxial[x0,radius,sigma,l,m,rinf,Mhalo,Rhalo,precision]";
	
wavepolarF::usage = 
		"Integrator for polar perturbations - full inhormogeneous integration wavepolarf[x0,radius,sigma,l,m,rinf,Mhalo,Rhalo,cs,cr,precision]";
		
waveaxialH::usage = 
		"Integrator for axial perturbations - hormogeneous  waveaxial[radius,l,m,rinf,Mhalo,Rhalo,precision]";


Begin["`Private`"];


{omega,drdrs,energy,angular} = Import["./src/geodesics.m"];

{mastera,sourcea,Xhora,dXhora,Xinfa,dXinfa}= Import["./src/final_axial.m"];

{masterav,sourceav,Xhorav,dXhorav,Xinfav,dXinfav} = Import["./src/final_axial_vac.m"];

{masterpf,sourcepf,H1horp,H0horp,Khorp,Whorp,\[Rho]horp,H1infp,dH1infp,H0infp,dH0infp,Kinfp,dKinfp,Winfp,dWinfp,\[Rho]infp,d\[Rho]infp} = Import["./src/final_polar_full.m"]

masteraxialH = Import["./src/final_axial_hom.m"];


(* ::Subsection:: *)
(*ODEs solver*)


(* ::Subsubsection::Closed:: *)
(*Axial inhomogeneous integrator*)


waveaxial[x0_,rp_,\[Sigma]\[Delta]_,l_,lm_,rinf_,Mhalo_,Rhalo_,precision_]:=Block[{\[Sigma]=\[Sigma]\[Delta],rpart=rp,
n,\[Mu],MBH=1,M=Mhalo,a0=Rhalo Mhalo,\[Epsilon]=10^-3,r,X,XH,dXH,XI,dXI,r1,r2,\[Omega],YH,rules,eqns,\[ScriptL]=l,\[ScriptM]=lm,Xh,X\[Infinity],ang,enr},

\[Mu]=1;

Xh[0]=x0;
X\[Infinity][0]=1;(* amplitude at infinity *)

enr=energy//.r->rp;
ang=angular//.r->rp;

PrintTemporary["\[Epsilon]=", ScientificForm[N[\[Epsilon],5]]];

rules={MaxSteps->Infinity,Method->"StiffnessSwitching",InterpolationOrder->All,WorkingPrecision->precision};

n=(\[ScriptL]-1)(\[ScriptL]+2)/2;

\[Omega]=\[ScriptM] omega//.r->rp;

{r1,r2}={2MBH(1+\[Epsilon]),Max[rinf/Abs[\[Omega]],2a0]};

eqns=If[Mhalo==0\[Or]Rhalo==0,masterav-sourceav==0,mastera-sourcea==0];

(* --- Boundary conditions at the horizon --- *)
{XH,dXH}=If[Mhalo==0\[Or]Rhalo==0,{Xhorav,dXhorav}//.r->r1,{Xhora,dXhora}//.r->r1];

(* --- Forward integration --- *)
YH=NDSolve[{eqns,X[r1]==XH,X'[r1]==dXH},X,{r,r1,r2},rules][[1]];

(* --- Boundary conditions at the infinity --- *)
{XI,dXI}=If[Mhalo==0\[Or]Rhalo==0,{Xinfav,dXinfav}//.r->r2,{Xinfa,dXinfa}//.r->r2];

PrintTemporary["|Jump| = ",Abs[((X[r2]//.YH)dXI-XI(X'[r2]//.YH))/XI]];

Return[{((X[r2]//.YH)dXI-XI(X'[r2]//.YH))/XI,1/(8\[Pi]) (\[ScriptL]+2)!/(\[ScriptL]-2)! Abs[X[r2]//.YH]^2}];

]


(* ::Subsubsection::Closed:: *)
(*Axial homogeneous integrator*)


waveaxialH[rp_,l_,lm_,rinf_,Mhalo_,Rhalo_,precision_]:=Block[{n,Xh0=1,X\[Infinity]0=1,\[Mu],
MBH=1,M=Mhalo,a0=Rhalo Mhalo,\[Epsilon]=10^-3,drs,r,X,XH,dXH,XI,dXI,r1,r2,\[Omega],YH,Y\[Infinity],rules,W,eqns,Xout,Xin,dEd\[Omega]g,\[ScriptL]=l,\[ScriptM]=lm,\[Delta],S\[Delta],Sd\[Delta],dEd\[Omega]gh,enr,ang},

\[Mu]=1;

n=(\[ScriptL]-1)(\[ScriptL]+2)/2;

drs=drdrs;

rules={MaxSteps->Infinity,Method->"StiffnessSwitching",InterpolationOrder->All,WorkingPrecision->precision};

PrintTemporary["(M,a0/M) = (",Mhalo,",",If[Mhalo==0,0,a0/Mhalo],")"];

enr=energy//.r->rp;
ang=angular//.r->rp;

\[Omega]=\[ScriptM] omega//.r->rp;

{r1,r2}={2MBH(1+\[Epsilon]),Max[rinf/Abs[\[Omega]],2a0]};

PrintTemporary["(r\[Infinity],r\[Infinity]/\[Omega],2a0) = (",N[r2 \[Omega]],",",N[rinf/\[Omega]],",",N[2a0],")"];

eqns=masteraxialH[[1]]==0;

(* ---- Boundary conditions at the horizon ---- *)

{XH,dXH}={masteraxialH[[3]],masteraxialH[[4]]}//.r->r1;

(* ---- Forward integration ---- *)

YH=NDSolve[{eqns,X[r1]==XH,X'[r1]==dXH},X,{r,r1,r2},rules][[1]];

(* ---- Boundary conditions at the infinity ---- *)

{XI,dXI}={masteraxialH[[5]],masteraxialH[[6]]}//.r->r2;

(* ---- Inward integration ---- *)

Y\[Infinity]=NDSolve[{eqns,X[r2]==XI,X'[r2]==dXI},X,{r,r2,r1},rules][[1]];

(* Define the worknskian *)

W=((X[r]//.YH)(X'[r]//.Y\[Infinity])-(X[r]//.Y\[Infinity])(X'[r]//.YH))/drs//.r->r2;

{S\[Delta],Sd\[Delta]}={Coefficient[masteraxialH[[2]],\[Delta][r]],Coefficient[masteraxialH[[2]],\[Delta]'[r]]};

{Xout,Xin}={W^-1 ((X[r]//.YH)S\[Delta] drs-D[(X[r]//.YH) drs Sd\[Delta],r]),W^-1 ((X[r]//.Y\[Infinity])S\[Delta] drs-D[(X[r]//.Y\[Infinity]) drs Sd\[Delta],r])}//.r->rp;

(* ---- Computing the fluxes ---- *)

{dEd\[Omega]g,dEd\[Omega]gh}=1/(8\[Pi]) (\[ScriptL]+2)!/(\[ScriptL]-2)! {Abs[Xout]^2,Abs[Xin]^2};


Return[{dEd\[Omega]g,dEd\[Omega]gh}]

]


(* ::Subsubsection:: *)
(*Polar inhomogeneous integrator*)


wavepolarF[x0_,rp_,\[Sigma]\[Delta]_,l_,lm_,rinf_,Mhalo_,Rhalo_,ccs_,ccsr_,precision_]:=Block[{cs=ccs,csr=ccsr,\[Sigma]=\[Sigma]\[Delta],rpart=rp,n,\[Mu],MBH=1,M=Mhalo,a0=Rhalo Mhalo,
\[Epsilon]=10^-4,r,r1,r2,\[Omega],YH,rules,eqns,\[ScriptL]=l,\[ScriptM]=lm,enr,ang,Kh,K\[Infinity]0,H0,H1,KK,W,\[Delta]\[Rho],H0H,KH,H1H,WH,\[Delta]\[Rho]H,KI,dKI},

Kh[0]=x0;
K\[Infinity]0=N[1,precision];(* amplitude at infinity *)

PrintTemporary["\[Epsilon]=", ScientificForm[N[\[Epsilon],5]]];

\[Mu]=1;

enr=energy//.r->rp;
ang=angular//.r->rp;

rules={MaxSteps->Infinity,Method->"StiffnessSwitching",InterpolationOrder->All,WorkingPrecision->precision};

n=(\[ScriptL]-1)(\[ScriptL]+2)/2;

\[Omega]=\[ScriptM] omega//.r->rp;

{r1,r2}={2MBH(1+\[Epsilon]),Max[rinf/Abs[\[Omega]],2a0]};

PrintTemporary[N[{rinf/Abs[\[Omega]],2a0},20]];

(* --- Boundary conditions at the horizon --- *)

{H1H,H0H,KH,WH,\[Delta]\[Rho]H}={H1horp,H0horp,Khorp,Whorp,\[Rho]horp}//.{r->r1};

eqns=Join[{masterpf[[1]]-\[Mu]*sourcepf[[1]]==0,masterpf[[2]]-\[Mu]*sourcepf[[2]]==0,masterpf[[3]]-\[Mu]*sourcepf[[3]]==0,
		masterpf[[4]]-\[Mu]*sourcepf[[4]]==0,masterpf[[5]]-\[Mu]*sourcepf[[5]]==0},{H1[r1]==H1H,H0[r1]==H0H,KK[r1]==KH,W[r1]==WH,\[Delta]\[Rho][r1]==\[Delta]\[Rho]H}];

(* --- Forward integration --- *)

YH=NDSolve[eqns,{H1,H0,KK,W,\[Delta]\[Rho]},{r,r1,r2},rules][[1]];

(* --- Boundary conditions at the infinity --- *)

{KI,dKI}={Kinfp,dKinfp}//.r->r2;

PrintTemporary["|jump| = ",Abs[((KK[r2]//.YH)dKI-KI(KK'[r2]//.YH))/KI]];

Return[{((KK[r2]//.YH)dKI-KI(KK'[r2]//.YH))/KI,1/(32\[Pi]) (\[ScriptL]+2)!/(\[ScriptL]-2)! Abs[KK[r2]//.YH]^2}];

]


(* ::Subsection:: *)
(*End*)


End[];

EndPackage[];
