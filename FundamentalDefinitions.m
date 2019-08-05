(* ::Package:: *)

BeginPackage["FundamentalDefinitions`"]

\[CapitalOmega]qqpp::usage="\[CapitalOmega]qqpp[ N] generates the symplectic form in basis (q1,q2,...,p1,p2,...) for N bosonic deg. of freedom"
\[CapitalOmega]qpqp::usage="\[CapitalOmega]qpqp[ N] generates the symplectic form in basis (q1,p1,q2,p2,...) for N bosonic deg. of freedom"

J::usage="J[ G] generates the complex structure associated with a cov. matrix G"
	
Transform::usage="Transform[ N] generates a N-dim. basis transformation matrix (q1,p1,q2,p2,...)\[Rule](q1,q2,...,p1,p2,...)"
	
SpBasis::usage="SpBasis[ N] generates sp(2N,R) Lie algebra basis"	
SpBasisNoUN::usage="spBasisNoUN[ N] generates sp(2N,R)/U(N) Lie algebra basis"	
SpBasisNoU1::usage="spBasisNoU1[ N] generates sp(2N,R)/U(1\!\(\*SuperscriptBox[\()\), \(N\)]\) Lie algebra basis"	
OBasis::usage="OBasis[ N] generates o(2N,R) Lie algebra basis"

MranSp::usage="MranSp[ N] generates a random NxN symplectic matrix"
MranO::usage="MranO[ N] generates a random NxN orthogonal matrix"

approxExp::usage="approxExp[ \[Epsilon], K]=(1 + \[Epsilon]/2 K)/(1 - \[Epsilon]/2 K) approximates MatrixExponential[\[Epsilon]X]"

logfunction::usage="logfunction[ x] takes value 0 when x=0 and log[x] otherwise"
ConditionalLog::usage"ConditionalLog[ x] calculates MatrixLog[ x] by applying logfunction to the eigenvalues of x"

invMSp::usage="invMSp[ m] inverts a symplectic matrix m as \!\(\*SuperscriptBox[\(m\), \(-1\)]\)=-\[CapitalOmega]m\[CapitalOmega]"

Begin["Private`"]

(* Symplectic forms *)
\[CapitalOmega]qqpp[n_]:=SparseArray[ArrayFlatten[{{0,IdentityMatrix[n]},{-IdentityMatrix[n],0}}]]; 
\[CapitalOmega]qpqp[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> ({{0, 1},{-1, 0}})]//SparseArray ;

(* Complex structure *)
J[G_]:=Module[{n,\[CapitalOmega]}, 
	n=Length[G]/2; \[CapitalOmega]=\[CapitalOmega]qpqp[n]; 
Return[-G.\[CapitalOmega]]
]

(* Basis transformation *)
Transform[n_]:=SparseArray[Join[Table[{k,2k-1}->1,{k,n}],Table[{n+k,2k}->1,{k,n}]]];

(* Lie algebra bases *)
SpBasis[n_]:=Module[{SYM,ASYM,notK,Tr,Trtran},
	SYM=Flatten[Table[SparseArray[{{i,j}->If[i!=j,N[1/Sqrt[2]],1],{j,i}->If[i!=j,N[1/Sqrt[2]],1]},{n,n}],{i,1,n},{j,i,n}],1];
	ASYM=Flatten[Table[SparseArray[{{i,j}->N[1/Sqrt[2]],{j,i}->-N[1/Sqrt[2]]},{n,n}],{i,1,n},{j,i+1,n}],1];
	notK=Join[
		Table[ArrayFlatten[{{m,0},{0,-m}}],{m,SYM}],
		Table[ArrayFlatten[{{m,0},{0,m}}],{m,ASYM}],
		Table[ArrayFlatten[{{0,m},{m,0}}],{m,SYM}],
		Table[ArrayFlatten[{{0,m},{-m,0}}],{m,SYM}]];
	Tr=Transform[n];Trtran=Transpose[Tr];
Return[Table[Trtran.Ki.Tr,{Ki,notK}]]
]
SpBasisNoUN[N_]:=Module[{Kold},
	Kold=SpBasis[N];
Return[Table[If[SymmetricMatrixQ[old],old],{old,Kold}]//DeleteCases[#,Null]&]
]
SpBasisNoU1[N_]:=Module[{Kold,notK1,notK2,Knew},
	Kold=SpBasis[N];
	notK1=Table[If[SymmetricMatrixQ[old]==False,old],{old,Kold}]//DeleteCases[#,Null]&;
	notK2=Table[If[Length[notK["NonzeroValues"]]==2,notK],{notK,notK1}]//DeleteCases[#,Null]&;
Return[DeleteCases[Kold,Alternatives@@notK2]]
]
OBasis[n_]:=Module[{notK,Tr,Trtran},
	notK=Flatten[Table[SparseArray[{{i,j}->1,{j,i}->-1},{2n,2n}],{i,1,2n},{j,i+1,2n}],1];
	Tr=Transform[n];Trtran=Transpose[Tr];
Return[Table[Trtran.Ki.Tr,{Ki,notK}]]
];

(* Matrix exp approximation *)
approxExp[t_,x_]:=Module[{dim,id},
	dim=Dimensions[x][[1]];id=IdentityMatrix[dim];
Return[(id+t/2 x).Inverse[(id-t/2 x)]]
]

(* Generate random Lie group elements *)
MranSp[n_]:=Module[{m1,m2,m3},
	m1=RandomReal[{0,1},{n,n}];
	m2=m1+Transpose[m1];
	m3=\[CapitalOmega]qpqp[n/2].m2;
Return[MatrixExp[m3]]
]
MranO[n_]:=Module[{m},
	m=RandomReal[{-1,1},{n,n}];
Return[Orthogonalize[m]]
]

(* Conditional logarithm *)
logfunction[x_]:=If[x==0,0,Log[x]]
ConditionalLog[x_]:=Module[{Diag,Tra,Mat},
	Diag=DiagonalMatrix[logfunction/@Eigenvalues[D.D]];
	Tra=Eigenvectors[D.D]//Transpose//Chop;
	Mat=Tra.Diag.Inverse[Tra]//Chop];

(* Inverting a symplectic matrix *)
invMSp[m_]:=Module[{\[CapitalOmega]}, \[CapitalOmega]=\[CapitalOmega]qpqp[Length[m]/2]; Return[-\[CapitalOmega].Transpose[m].\[CapitalOmega]]]

End[]

EndPackage[]



