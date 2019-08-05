(* ::Package:: *)

BeginPackage["Applications`",{"FundamentalDefinitions`"}]

CoPBos::usage="Generates function which calculates complexity w.r.t. target state JT"
CoPgradBos::usage="Generates function which calculates complexity gradient w.r.t. target state JT on basis K"

RestrictionAA::usage="Returns the resAA parameter for a given partition of a system"
EoPBos::usage="Generates function which calculates EoP w.r.t. given restriction (for bosons)"
EoPgradBos::usage="Generates function which calculates EoP gradient w.r.t. given restriction (for bosons)"
EoPFerm::usage="Generates function which calculates EoP w.r.t. given restriction (for fermions)"
EoPgradFerm::usage="Generates function which calculates EoP gradient w.r.t. given restriction (for fermions)"

energyBos::usage="Generates energy function for bosonic quadr. Hamiltonian with h"
energygradBos::usage="Generates energy gradient function for bosonic quadr. Hamiltonian with h"

energyFerm::usage="Generates energy function for fermionic quadr. Hamiltonian with h"
energygradFerm::usage="Generates energy gradient function for fermionic quadr. Hamiltonian with h"

Begin["Private`"]

(* Complexity of purification - function *)
CoPBos[JT_]:=Function[{M,Jref0},Module[{invM,D},
	D=M.Jref0.invMSp[M].Inverse[JT]; 
	Sqrt[Total[Log[#1]^2&/@Eigenvalues[D]]/8]]
];

(* Complexity of purification - gradient *)
CoPgradBos[JT_]:=Function[{M,Jref0,K},Module[{dimA,dimB,invM,invJT,D,invD,dD},
	dimB=Length[K[[1]]]; dimA=Length[M]-dimB;
	invM=invMSp[M]; invJT=Inverse[JT];
	D=M.Jref0.invM.Inverse[JT]; invD=Inverse[D];
	dD=Table[M.(KIi.Jref0-Jref0.KIi).invM.invJT,{KIi,K}];
	Table[2Tr[MatrixLog[D].invD.dDi],{dDi,dD}]]
];

(* Function to generate restriction - for use in EoP *)
RestrictionAA[partition_]:=Module[{n,n1,n2,A1,B1,A2,B2},
	{A1,B1,A2,B2}=partition;
	n=A1+B1+A2+B2; n1=A1+B1;
	Join[Table[x,{x,2A1}],Table[y+2n1,{y,2A2}]]
];

(* Entanglement of purification - function *)
EoPBos[ResAA_]:=Function[{M,J0}, Module[{n,D,logDD},
	n=Length[J0];
	D=(M.(IdentityMatrix[n]+I J0).invMSp[M]/2)[[ResAA,ResAA]];
	logDD=ConditionalLog[D.D];
	Re[1/2 Tr[D.logDD]]//Chop]
];
EoPFerm[ResAA_]:=Function[{M,J0}, Module[{n,D,logD},
	n=Length[J0];
	D=(M.(IdentityMatrix[n]+I J0).invMSp[M]/2)[[ResAA,ResAA]];
	logD=ConditionalLog[D];
	Re[-Tr[D.Mat]]//Chop]
];

(* Entanglement of purification - gradient *)
EoPgradBos[ResAA_]:=Function[{M,J0,K},Module[{n,nn,invM,dJ,D,dD,logDD},
	n=Length[J0];
	invM=invMSp[M];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,K}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logDD=ConditionalLog[D.D];
	Table[Re[1/2 Tr[dDi.logDD]]//Chop,{dDi,dD}]]
];
EoPgradFerm[ResAA_]:=Function[{M,J0,K},Module[{n,nn,invM,KI,dJ,D,dD,logD},
	n=Length[J0]; nn=Length[J0];
	invM=invMSp[M];
	KI=Table[PadLeft[Ki,{n,n}],{Ki,K}];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,KI}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logD=ConditionalLog[D];
	Table[Re[-Tr[dDi.logD]]//Chop,{dDi,dD}]]
];

(* Energy for bososnic d.o.f. - function *)
energyBos[h_]:=Function[{M,J0},Module[{\[CapitalOmega]0,G},
	\[CapitalOmega]0=\[CapitalOmega]qpqp[Length[M]/2];
	G=M.J0.\[CapitalOmega]0.Transpose[M]; 1/4 Tr[h.G]]
];

(* Energy for bososnic d.o.f. - gradient *)
energygradBos[h_]:=Function[{M,J0,K},Module[{\[CapitalOmega]0},
	\[CapitalOmega]0=\[CapitalOmega]qpqp[Length[M]/2];
	Table[1/4 Tr[h.M.(KK.J0.\[CapitalOmega]0+J0.\[CapitalOmega]0.Transpose[KK]).Transpose[M]],{KK,K}]]
];

(* Energy for fermionic d.o.f. - function *)
energyFerm[h_]:=Function[{M,J0},Module[{J},
	J=M.J0.Inverse[M]; 1/4 Tr[h.J]]
];

(* Energy for fermionic d.o.f. - gradient *)
energygradFerm[h_]:=Function[{M,J0,K},
Table[1/4 Tr[h.M.(KK.J0-J0.KK).Inverse[M]],{KK,K}]
];

End[]

EndPackage[]



