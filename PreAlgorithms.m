(* ::Package:: *)

BeginPackage["PreAlgorithms`",{"FundamentalDefinitions`"}]

ExtractStdForm::usage="ExtractStdForm[ G] 
	\!\(\*
StyleBox[\"Input\",\nFontWeight->\"Bold\"]\): Arbitrary covariance matrix G (bosonic/fermionic) 
	\!\(\*
StyleBox[\"Output\",\nFontWeight->\"Bold\"]\): ri parameters required to construct the std. form, as well as the transformation that brought CM into this form"
 
ConstructPurificationBos::usage="ConstructPurificationBos[ r, {\!\(\*SubscriptBox[\(N\), \(A\)]\),\!\(\*SubscriptBox[\(N\), \(B\)]\)}] 
	\!\(\*
StyleBox[\"Input\",\nFontWeight->\"Bold\"]\): List of ri parameters, dimensions of subsystems A and B
	\!\(\*
StyleBox[\"Output\",\nFontWeight->\"Bold\"]\): Complex structure associated with purified state"

ConstructPurificationFerm::usage="ConstructPurificationFerm[ r, {\!\(\*SubscriptBox[\(N\), \(A\)]\),\!\(\*SubscriptBox[\(N\), \(B\)]\)}] 
	\!\(\*
StyleBox[\"Input\",\nFontWeight->\"Bold\"]\): List of ri parameters, dimensions of subsystems A and B
	\!\(\*
StyleBox[\"Output\",\nFontWeight->\"Bold\"]\): Complex structure associated with purified state"

Begin["Private`"]

ExtractStdForm[G_]:=Module[{NT=Dimensions[G][[1]]/2,\[CapitalOmega]T,J,RI,SWT,tra,resc,Tra,GAB,ON,TRA,check,checklist,DiagElems,Diag,Mtra,rlist}, 
	\[CapitalOmega]T=\[CapitalOmega]qpqp[NT]; J=-G.\[CapitalOmega]T;

	RI=SparseArray[{{x_,x_}->-I Abs[Re[I^x]]+ Abs[Re[I^(x+1)]],{x_,y_}/;y-x==1->Abs[Re[I^y]],{x_,y_}/;x-y==1->I Abs[Re[I^x]]},{2NT,2NT}];
	SWT=SparseArray[{{x_,y_}/;y-x==1->Abs[Re[I^y]],{x_,y_}/;x-y==1-> Abs[Re[I^x]]},{2NT,2NT}]; 
	tra=Transpose[Eigenvectors[J]].Transpose[RI]//Chop;  
	resc=SWT.Transpose[tra].Inverse[\[CapitalOmega]T].tra//Chop//Abs//Sqrt//Inverse//Chop; 
	Tra=Inverse[tra.resc]//Chop ; 
	GAB=(Tra.G.Transpose[Tra])//Chop; 
	ON=Eigenvectors[GAB];
	TRA=ON.Tra;

	check=TRA.\[CapitalOmega]T.Transpose[TRA]//Chop;
	checklist=(check.ConstantArray[1,2 NT])[[;;;;2]]//Round;
	DiagElems=Replace[checklist,Dispatch[{-1->a,1->b}],1];
	Diag=ArrayFlatten[DiagonalMatrix[DiagElems]/.{a->{{0,1},{1,0}},b->IdentityMatrix[2]}];
	Mtra=Diag.TRA;

	rlist=Diagonal[1/2 ArcCosh[((TRA.G.Transpose[TRA])//Chop)]][[1;;-1;;2]];
Return[{rlist,Mtra}]
]

(* Purified state for bosons *)
ConstructPurificationBos[rlist_,{n_,m_}]:=Module[{cosh,sinh,diag,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0},
	cosh=Flatten[Table[{Cosh[2 rr],Cosh[2 rr]},{rr,rlist}]];
	sinh=Flatten[Table[{Sinh[2 rr],-Sinh[2 rr]},{rr,rlist}]];

	Q14=DiagonalMatrix[Hold/@cosh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sinh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q5=If[m==n,Null,IdentityMatrix[2(m-n)]//SparseArray];

	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{Q23,Q14}}]//SparseArray, G0=ArrayFlatten[{{Q14,Q23,0},{Q23,Q14,0},{0,0,Q5}}]//SparseArray];

	\[CapitalOmega]0=\[CapitalOmega]qpqp[m+n];
Return[-G0.\[CapitalOmega]0]
];

(* Purified state for fermions *)
ConstructPurificationFerm[rlist_,{n_,m_}]:=Module[{cos,sin,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0},
	cos=Table[ArrayFlatten[{{0,Cos[2 rr]},{-Cos[2 rr],0}}],{rr,rlist}];
	sin=Table[ArrayFlatten[{{0,Sin[2 rr]},{Sin[2 rr],0}}],{rr,rlist}];

	Q14=DiagonalMatrix[Hold/@cos]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sin]//ReleaseHold//ArrayFlatten//SparseArray;
	
	Q5=If[m==n,Null,ArrayFlatten[DiagonalMatrix[ConstantArray[a,m-n]]/.a -> {{0, 1},{-1, 0}}]];
	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{-Q23,Q14}}]//SparseArray, 
	G0=ArrayFlatten[{{Q14,Q23,0},{-Q23,Q14,0},{0,0,Q5}}]//SparseArray];
	
	\[CapitalOmega]0=\[CapitalOmega]qpqp[m+n];
Return[-G0.\[CapitalOmega]0]
]

End[]

EndPackage[]






