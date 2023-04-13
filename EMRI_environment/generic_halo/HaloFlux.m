(* ::Package:: *)

BeginPackage["HaloFlux`"];


(* ::Subsubsection:: *)
(*Usage*)


background::usage = 
		"Computes the background metric functions for a given mass profile";

As\[Psi]h::usage = 
		"Computes the asymptotic condition at horizon";
		
As\[Psi]\[Infinity]::usage = 
		"Computes the asymptotic condition at infinity";

\[Rho]profile::usage = 
		"Provides the density profile";

FluxGen::usage = 
		"Provides the density profile";


Begin["`Private`"];


{soldm,soldF,rulesG}=Import["./src/rulesmetric.m"];
{rulesGeo,rulesEL}=Import["./src/rulesgeodesics.m"];
{master\[Psi],source\[Psi],master\[Psi]log,source\[Psi]log,master\[Psi]H,master\[Psi]Hlog}=Import["./src/master_axial.m"];


year = 365;
day = 3600 24;
GN = 6674 10^-14;
cc = 299792458 10^-3;
msun = 14768 10^-4;
Rs = 1495978707/10;
pc = 96939420213600/\[Pi];
\[Gamma]E = 57721 10^-5;


(* ::Subsubsection::Closed:: *)
(*Halo density profile*)


\[Rho]profile[OptionsPattern[{prec\[Rho]->25}],r_,Mhalo_,param_,mode_]:=Block[{rh=2,\[Rho],a0,rc,rules,fact,xx},

rules = {WorkingPrecision->OptionValue[prec\[Rho]],MaxRecursion->1000};
fact[xx_] = (1-rh/xx);

a0=If[ArrayQ[param],param[[1]],param];
rc=If[ArrayQ[param],param[[2]],0];

(*--- rc is the cut-off for NFW, while the \[Infinity] for Einasto ---*)
(*--- for Einasto a0 = re ---*)

If[mode!=("HQ"\[Or]"NFW"\[Or]"Ein"),Print["wrong profile choice"];Return[],Continue[]];

\[Rho]=Which[

	mode=="HQ",Mhalo (a0+2)/(2\[Pi] r (r+a0)^3) fact[r],

	mode=="NFW",

	Mhalo/NIntegrate[4\[Pi] xx^2  a0/xx 1/(1+xx/a0)^2 fact[xx],{xx,rh,rc},WorkingPrecision->OptionValue[prec\[Rho]],MaxRecursion->1000] a0^3/(r (r+a0)^2) fact[r] UnitStep[rc-r],

	mode=="Ein",

	Mhalo/NIntegrate[4\[Pi] xx^2  Exp[(-53/3)((xx/a0)^(1/6)-1)]fact[xx],{xx,rh,rc},WorkingPrecision->OptionValue[prec\[Rho]],MaxRecursion->1000] Exp[(-53/3)((r/a0)^(1/6)-1)]fact[r],

	mode=="HQB", (-((2.9732001769144445`*^-18 (1-8/r)^1.5001022653136356`)/r^1.5000599128406606`)+3913/(5000000000000000000000 (1+r/910000000000)^3 r)) HeavisideTheta[-223683.53464147683823920250163736371568867415`15.+r]+(2.9732001769144445`*^-18 (1-8/r)^1.5001022653136356` HeavisideTheta[-8+r])/r^1.5000599128406606`,

	mode=="HQcust",N[(-((2.9550986187354686`*^-12 (1-8/r)^1.5000477723165158`)/r^1.500431626921895`)+3913/(50000000000000 (1+r/91000000)^3 r)) HeavisideTheta[-22.36835429299709075401730962756152630075`15.+r]+(2.9550986187354686`*^-12 (1-8/r)^1.5000477723165158` HeavisideTheta[-8+r])/r^1.500431626921895`,\[Infinity]]
];

Return[N[\[Rho],OptionValue[prec\[Rho]]]]

]


(* ::Subsubsection::Closed:: *)
(*Background geometry integrator*)


background[OptionsPattern[{precision->25,eps->10^-4}],mode_,Mhalo_,param_,rout_]:=Block[{rh,r1,rules,m,F,\[Rho],r,
M=1,Nsolm,NsolF,a0,drt,rt,Nsolrt,solrISCO,solrLR,\[CapitalOmega]ISCO,\[CapitalOmega]LR,u\[Phi],ut,BCdrt,rinf=rout,eqns,BCs,solmetric,extraprec},

extraprec = 40;

rh=2M;
r1=2M(1+OptionValue[eps]);

rules={MaxSteps->10^6,Method -> "StiffnessSwitching",InterpolationOrder->All};

a0=If[ArrayQ[param],param[[1]],param];

\[Rho][r_]=\[Rho]profile[prec\[Rho]->OptionValue[precision]+extraprec,r,Mhalo,param,mode];

(*--- Integrate the mass function m[r] ---*)

Nsolm=NDSolve[{m'[r]==4 \[Pi] r^2 \[Rho][r],m[rh]==M},m,{r,rh,rinf},Join[rules,{WorkingPrecision->OptionValue[precision]+extraprec}]][[1]];

(*--- Integrate F[r] and rt[r] ---*)

eqns = {F'[r]==(2 F[r] m[r])/(r-2 m[r])/r,rt'[r]==1/Sqrt[F[r](1-2m[r]/r)]}//.Nsolm;

BCs = {F[rinf]==1-2m[rinf]/rinf//.Nsolm,rt[rinf]==rinf+2 m[rinf] Log[rinf/2(m[rinf])-1]//.Nsolm};

NsolF=NDSolve[Join[eqns,BCs],{F,rt},{r,rinf,r1},Join[rules,{WorkingPrecision->OptionValue[precision]+extraprec/2}]][[1]];

solmetric=Join[NsolF,Nsolm];

(*--- geodetic properties ---*)

solrISCO=FindRoot[(r m[r]-6m[r]^2+r^2 m'[r]==0)//.solmetric,{r,6},WorkingPrecision->OptionValue[precision]][[1,2]];

solrLR=FindRoot[r==3m[r]/.solmetric,{r,3},WorkingPrecision->OptionValue[precision]][[1,2]];

{\[CapitalOmega]LR,\[CapitalOmega]ISCO}={Sqrt[F[solrLR]]/solrLR,Sqrt[(F[solrISCO] m[solrISCO])/(solrISCO^2-2 solrISCO m[solrISCO])] 1/Sqrt[solrISCO]}//.solmetric;

Return[{solmetric,{solrLR,solrISCO,\[CapitalOmega]LR,\[CapitalOmega]ISCO}}]

]


(* ::Subsubsection::Closed:: *)
(*Asymptotic expansion at the horizon*)


As\[Psi]h[OptionsPattern[{precision->25}],order_,l_,\[CapitalOmega]_,rr_,solmetric_]:=
Block[{rh,M=1,\[ScriptL]=l,\[Omega]=\[CapitalOmega],\[Beta],m,F,mcoef,Fcoef,\[Psi],r,rulesmcoef,rulesFcoef,seriesh,drt,rulesRt,rt,d\[Psi]h},

(*--- Tortoise coordinates ---*)

rh=2M;

drt=1/Sqrt[F[r](1-2m[r]/r)];

rulesRt={rt[r]->0,rt'[r]->drt,rt''[r]->D[drt,r]};

(*--- Getting m and F expansion coef with their numeric interpolating functions --- *)

rulesmcoef=Table[mcoef[i]->D[m[r],{r,i}]/Factorial[i]//.solmetric,{i,2,order}];
rulesFcoef=Table[Fcoef[i]->D[F[r],{r,i}]/Factorial[i]//.solmetric,{i,1,order}];

(*--- Expansion of the function at the horizon ---*)

m[r_]=M+Sum[mcoef[i](r-rh)^i,{i,2,order}];
F[r_]=Sum[Fcoef[i](r-rh)^i,{i,1,order}];
\[Psi][r_]=Sum[Exp[-I \[Omega] rt[r]]\[Beta][i](r-rh)^i,{i,0,order}];

(*--- Expansion of the master equation ---*)

\[Beta][0]=1;

seriesh=Assuming[r>rh,Collect[Series[master\[Psi]//.rulesRt,{r,rh,order}],r-rh]];

(*--- Getting all the coefficients ---*)
Table[\[Beta][i]=Solve[(Coefficient[seriesh,r-rh,i]==0)//.Join[rulesmcoef,rulesFcoef],\[Beta][i],WorkingPrecision->OptionValue[precision]][[1,1,2]]//.r->rr,{i,1,order}];

d\[Psi]h=\[Psi]'[r]//.Join[rulesmcoef,rulesFcoef,rulesRt];

PrintTemporary["+Asymptotic horizon computed"];

Return[{\[Psi][rr]//.rt[rr]->0,d\[Psi]h//.r->rr}];
Return[{\[Psi][rr]//.solmetric,Exp[-I \[Omega] rt[rr]]d\[Psi]h//.solmetric//.r->rr}];

]


(* ::Subsubsection::Closed:: *)
(*Asymptotic expansion at infinity*)


As\[Psi]\[Infinity][OptionsPattern[{precision->25}],order_,l_,\[CapitalOmega]_,rr_,solmetric_]:=
Block[{M=1,\[ScriptL]=l,\[Omega]=\[CapitalOmega],\[Gamma],m,F,mcoef,Fcoef,\[Psi],r,series\[Infinity],hyp,drt,ddrt,rulesRt,rt,d\[Psi]\[Infinity],rulesmcoef,rulesFcoef},

(*--- Tortoise coordinates ---*)

drt=1/Sqrt[F[r](1-2m[r]/r)];

rulesRt={rt[r]->0,rt'[r]->drt,rt''[r]->D[drt,r]};

(*--- At \[Infinity] we assume everything is as Schwarzschild with M+Mhalo ---*)

rulesmcoef=Table[mcoef[i]->D[m[r],{r,i}]/Factorial[i]//.solmetric,{i,0,order}];
rulesFcoef=Table[Fcoef[i]->D[F[r],{r,i}]/Factorial[i]//.solmetric,{i,0,order}];

m[r_]=Sum[mcoef[i]r^(-i),{i,0,order}];
F[r_]=Sum[Fcoef[i]r^(-i),{i,0,order}];

\[Psi][r_]=Sum[Exp[I \[Omega] rt[r]]\[Gamma][i]r^(-i),{i,0,order}];

\[Gamma][0]=1;

series\[Infinity]=Normal[Series[master\[Psi]//.rulesRt,{r,\[Infinity],order+2}]];

(*--- Getting all the coefficients ---*)
Table[\[Gamma][i]=Solve[(Coefficient[series\[Infinity],r,-i-1]==0)//.Join[rulesmcoef,rulesFcoef],\[Gamma][i],WorkingPrecision->OptionValue[precision]][[1,1,2]]//.r->rr,{i,1,order}];

d\[Psi]\[Infinity]=\[Psi]'[r]//.Join[rulesmcoef,rulesFcoef,rulesRt];

PrintTemporary["+Asymptotic infinity computed"];

Return[{\[Psi][rr]//.rt[rr]->0,d\[Psi]\[Infinity]//.r->rr}];

Return[{\[Psi][rr]//.solmetric,Exp[I \[Omega] rt[rr]]d\[Psi]\[Infinity]//.solmetric//.r->rr}];


]


(* ::Subsubsection:: *)
(*Axial flux generator*)


FluxGen[OptionsPattern[{prec->30,epsilon->10^-4,order->4,rinf->10^13}],l_,lm_,rp_,robs_,mode_,Mhalo_,param_]:=
Block[{\[Mu]=1,M=1,\[Omega],\[ScriptL]=l,\[ScriptM]=lm,r2,rh,a0,\[Psi],r,r1,rules,\[Rho],ODE\[Psi],ODE\[Psi]t,\[Psi]rh,d\[Psi]rh,sol\[Psi]H,\[Psi]r\[Infinity],d\[Psi]r\[Infinity],sol\[Psi]\[Infinity],
coef,coef\[Delta],coefd\[Delta],dcoefd\[Delta],solg,m,F,rt,drt,Wrt,source,\[Psi]out,\[Theta],\[Phi],\[Delta],r\[Infinity]=OptionValue[rinf],\[Psi]\[Infinity],\[Psi]H,cutoff},

rh=2M;
a0=If[ArrayQ[param],param[[1]],param];

cutoff=If[mode=="NFW",param[[2]],0];

rules={MaxSteps->10^6,InterpolationOrder->All,WorkingPrecision->OptionValue[prec],Method->"StiffnessSwitching"};

(*--- Call the routine for the background solutions ---*)

solg=background[precision->OptionValue[prec],eps->OptionValue[epsilon],mode,Mhalo,param,r\[Infinity]][[1]];

PrintTemporary["Metric precision = ",Precision[solg]];

drt=1/Sqrt[F[r](1-2m[r]/r)]//.solg;
\[Omega]=lm*Sqrt[(F[rp] m[rp])/(rp^2-2 rp m[rp])] 1/Sqrt[rp]//.solg;

{r1,r2}={rh(1+OptionValue[epsilon]),Max[robs/Abs[\[Omega]],2Max[a0,cutoff]]};
Print[ScientificForm[r2]];
If[r2>r\[Infinity],Print["domain of integration not large enough"]; Return[]];

ODE\[Psi]=(master\[Psi]==0)//.solg;

PrintTemporary["Integration at the horizon"];

(*--- Forward integration ---*)
{\[Psi]rh,d\[Psi]rh}=As\[Psi]h[precision->Precision[{\[Omega],solg}],OptionValue[order],l,\[Omega],r1,solg];

sol\[Psi]H=NDSolve[{ODE\[Psi],\[Psi][r1]==\[Psi]rh,\[Psi]'[r1]==d\[Psi]rh},\[Psi],{r,r1,r2},rules][[1]];

\[Psi]H[r_]=\[Psi][r]//.Join[sol\[Psi]H,solg];

PrintTemporary["Integration at infinity"];

(*--- Backckward integration---*)
{\[Psi]r\[Infinity],d\[Psi]r\[Infinity]}=As\[Psi]\[Infinity][precision->Precision[{\[Omega],solg}],OptionValue[order],l,\[Omega],r2,solg];

(*d\[Psi]r\[Infinity] = d\[Psi]r\[Infinity]-0I*\[Omega]*drt*\[Psi]r\[Infinity]//.r->r2//.solg;
ODE\[Psi]t=(master\[Psi]==0)//.solg;*)

sol\[Psi]\[Infinity]=NDSolve[{ODE\[Psi],\[Psi][r2]==\[Psi]r\[Infinity],\[Psi]'[r2]==d\[Psi]r\[Infinity]},\[Psi],{r,r2,r1},rules][[1]];
\[Psi]\[Infinity][r_]=\[Psi][r]//.Join[sol\[Psi]\[Infinity],solg];

(*--- Inomogenehous integration---*)

Wrt=(\[Psi]H[r]*\[Psi]\[Infinity]'[r]-\[Psi]H'[r]*\[Psi]\[Infinity][r])/drt//.r->rp;

coef=\[Psi]H[r] drt//.sol\[Psi]H;

{coef\[Delta],dcoefd\[Delta]}={Coefficient[source\[Psi],\[Delta][r]],Coefficient[source\[Psi],\[Delta]'[r]]};

\[Psi]out=\[Psi]\[Infinity][r2](coef*coef\[Delta]/Wrt-D[coef dcoefd\[Delta],r]/Wrt)//.Join[{r->rp,\[Theta]->\[Pi]/2,\[Phi]->0},solg,rulesEL];

Return[1/(8\[Pi]) (l+2)!/(l-2)! Abs[\[Psi]out]^2]

]


(* ::Subsection::Closed:: *)
(*End*)


End[];

EndPackage[];
