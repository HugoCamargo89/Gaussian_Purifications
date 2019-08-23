(* ::Package:: *)

BeginPackage["GaussianOptimization`"];


$GOFunctions=Sort[{"GO\[CapitalOmega]qqpp","GO\[CapitalOmega]qpqp","GOTransformtoJ","GOTransform","GOMranSp","GOMranO","GOSpBasis","GOSpBasisNoUN","GOSpBasisNoU1","GOOBasis","GOapproxExp","GOlogfunction","GOConditionalLog","GOinvMSp","GOOptimize","GOExtractStdForm","GOConstructPurificationBos","GOConstructPurificationFerm"}];


$GOApplicationsFunctions=Sort[{"GORestrictionAA","GOEoPBos","GOEoPFerm","GOEoPgradBos","GOEoPgradFerm","GOCoPBos","GOCoPgradBos"}];


(* Optimization algorithm *)
GOOptimize::usage="Optimization of a scalar function over the manifold of Gaussian states. 

The algorithm initializes at at various specified starting points in the manifold, each of which defines a trajectory for the optimization. 
At each multiple of 5 steps, the algorithm retains only the 20% of the trajectories with the lowest function values and highest gradient norm, 
respectively. The optimization is terminated when the final remaining trajectory (or trajectories) first satisfies one of three stopping criteria: 
a lower bound on the gradient norm, a limit on the number of iterations, and a relative difference in function value below \!\(\*SuperscriptBox[\(10\), \(-10\)]\). 

\!\(\*
StyleBox[\"Input\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"arguments\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"are\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) function to be optimised f[ M, J0] 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"2\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) associated gradient function df[ M, J0, Kbasis] 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"3\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) initial J 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"4\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) list of initial transformations {M1,M2,...} 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"5\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) specification of basis types for subsystems A and B { 'basisA', 'basisB'} (must be strings and must specify 'None' if no optimization occurs in one subsystem) 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"6\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) no. of degrees of freedom in subsystems { dimA, dimB} (must specify dimB=0 and basisB='None' if there is no splitting into subsystems) 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"7\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) tolerance on norm of gradient as stopping criterion
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"8\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) iteration limit as stopping criterion (default=infinity)

\!\(\*
StyleBox[\"Output\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"arguments\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"are\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\":\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\" \",\nFontSlant->\"Italic\"]\)
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) f\!\(\*
StyleBox[\"inal\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"function\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"value\",\nFontWeight->\"Plain\"]\)
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"2\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) final transformation M 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"3\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) list of the no. of corrections of the step size per iteration 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"4\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) list of the function values for the final surviving trajectory 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"5\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) list of the values of the gradient norm for the final surviving trajectory 
\!\(\*
StyleBox[\"(\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"6\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\")\",\nFontWeight->\"Bold\"]\) \!\(\*
StyleBox[\"Natural\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"metric\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"on\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"the\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\"manifold\",\nFontWeight->\"Plain\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Plain\"]\)
";


(* Fundamental definitions and basic tools *)
(* Symplectic forms and various matrices *)
GO\[CapitalOmega]qqpp::usage="GO\[CapitalOmega]qqpp[ N] generates the symplectic form in basis (q1,q2,...,p1,p2,...) for N bosonic deg. of freedom";


GO\[CapitalOmega]qpqp::usage="GO\[CapitalOmega]qpqp[ N] generates the symplectic form in basis (q1,p1,q2,p2,...) for N bosonic deg. of freedom";


GOTransformtoJ::usage="GOTransformtoJ[ G] generates the complex structure associated with a cov. matrix G by computing J=-G.\[CapitalOmega]";


GOTransform::usage="GOTransform[ N] generates a N-dim. basis transformation matrix (q1,p1,q2,p2,...)\[Rule](q1,q2,...,p1,p2,...)";


GOMranSp::usage="MranSp[ N] generates a random NxN symplectic matrix";


GOMranO::usage="MranO[ N] generates a random NxN orthogonal matrix";


(* Lie algebra bases *)
GOSpBasis::usage="SpBasis[ N] generates sp(2N,R) Lie algebra basis";


GOSpBasisNoUN::usage="spBasisNoUN[ N] generates sp(2N,R)/U(N) Lie algebra basis";


GOSpBasisNoU1::usage="spBasisNoU1[ N] generates sp(2N,R)/U(1\!\(\*SuperscriptBox[\()\), \(N\)]\) Lie algebra basis";


GOOBasis::usage="OBasis[ N] generates o(2N,R) Lie algebra basis";


(* Purifications and the standard form *)
GOExtractStdForm::usage="GOExtractStdForm[ G] returns the parameters for constructing the standard form of G (list of the \!\(\*SubscriptBox[\(r\), \(i\)]\) and the transformation matrix that puts G into standard form)";


GOConstructPurificationBos::usage="GOConstructPurificationBos[ rlist, { dimA, dimB}] generates the complex structure for the (dimA+dimB)-dimensional purification of a dimA-dimensional bosonic Gaussian state in standard form";


GOConstructPurificationFerm::usage="GOConstructPurificationBos[ rlist, { dimA, dimB}] generates the complex structure for the (dimA+dimB)-dimensional purification of a dimA-dimensional fermionic Gaussian state in standard form";


(* Tools for computational efficiency *)
GOapproxExp::usage="approxExp[ \[Epsilon], K]=(1 + \[Epsilon]/2 K)/(1 - \[Epsilon]/2 K) approximates MatrixExponential[\[Epsilon]X]";


GOlogfunction::usage="logfunction[ x] takes value 0 when x=0 and log[x] otherwise";


GOConditionalLog::usage="GOConditionalLog[ x] calculates MatrixLog[ x] by applying logfunction to the eigenvalues of x";


GOinvMSp::usage="invMSp[ m] inverts a symplectic matrix m as \!\(\*SuperscriptBox[\(m\), \(-1\)]\)=-\[CapitalOmega]m\[CapitalOmega]";


(* Application: Entanglement of purification *)
GORestrictionAA::usage="GORestrictionAA[ {dimA1, dimB1, dimA2, dimB2}] generates the restriction function for use in calculating EoP for a system decomposed into degrees of freedom {dimA1, dimB1, dimA2, dimB2}";


GOEoPBos::usage="GOEoPBos[ Restriction] generates a function f[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\)] which calculates the bosonic EoP for the given restriction";


GOEoPFerm::usage="GOEoPBos[ Restriction] generates a function f[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\)] which calculates the fermionic EoP for the given restriction";


GOEoPgradBos::usage="GOEoPgradBos[ Restriction] generates a function df[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\), K] which calculates the gradient of the bosonic EoP for the given restriction with respect to the Lie algebra basis K";


GOEoPgradFerm::usage="GOEoPgradFerm[ Restriction] generates a function df[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\), K] which calculates the gradient of the fermionic EoP for the given restriction with respect to the Lie algebra basis K";


(* Application: Complexity of purification *)
GOCoPBos::usage="GOCoPBos[ \!\(\*SubscriptBox[\(J\), \(T\)]\)] generates a function f[ M, \!\(\*SubscriptBox[\(J\), \(R0\)]\)] which calculates the bosonic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\)";


GOCoPgradBos::usage="GOCoPBos[ \!\(\*SubscriptBox[\(J\), \(T\)]\)] generates a function df[ M, \!\(\*SubscriptBox[\(J\), \(R0\)]\), K] which calculates the gradient of the bosonic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\) and the Lie algebra basis K";


GOenergyBos::usage="GOenergyBos[ h] generates an energy function f[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\)] for a bosonic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\)";


GOenergygradBos::usage="GOenergygradBos[ h] generates an energy gradient function df[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\), K] for a bosonic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\) with respect to the Lie algebra basis K";


GOenergyFerm::usage="GOenergyFerm[ h] generates an energy function f[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\)] for a fermionic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\)";


GOenergygradFerm::usage="GOenergygradBos[ h] generates an energy gradient function df[ M, \!\(\*SubscriptBox[\(J\), \(0\)]\), K] for a fermionic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\) with respect to the Lie algebra basis K";


Begin["Private`"];

(* -----------------------------------------------------------Fundamental definitions and basic tools--------------------------------------------------------------- *)
(* Symplectic forms *)
GO\[CapitalOmega]qqpp[n_]:=SparseArray[ArrayFlatten[{{0,IdentityMatrix[n]},{-IdentityMatrix[n],0}}]]; 
GO\[CapitalOmega]qpqp[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> ({{0, 1},{-1, 0}})]//SparseArray ;

(* Complex structure *)
GOTransformtoJ[G_]:=Module[{n,\[CapitalOmega]}, 
	n=Length[G]/2; \[CapitalOmega]=GO\[CapitalOmega]qpqp[n]; 
Return[-G.\[CapitalOmega]]];

(* Basis transformation *)
GOTransform[n_]:=SparseArray[Join[Table[{k,2k-1}->1,{k,n}],Table[{n+k,2k}->1,{k,n}]]];

(* Lie algebra bases *)
GOSpBasis[n_]:=Module[{SYM,ASYM,notK,Tr,Trtran},
	SYM=Flatten[Table[SparseArray[{{i,j}->If[i!=j,N[1/Sqrt[2]],1],{j,i}->If[i!=j,N[1/Sqrt[2]],1]},{n,n}],{i,1,n},{j,i,n}],1];
	ASYM=Flatten[Table[SparseArray[{{i,j}->N[1/Sqrt[2]],{j,i}->-N[1/Sqrt[2]]},{n,n}],{i,1,n},{j,i+1,n}],1];
	notK=Join[
		Table[ArrayFlatten[{{m,0},{0,-m}}],{m,SYM}],
		Table[ArrayFlatten[{{m,0},{0,m}}],{m,ASYM}],
		Table[ArrayFlatten[{{0,m},{m,0}}],{m,SYM}],
		Table[ArrayFlatten[{{0,m},{-m,0}}],{m,SYM}]];
	Tr=GOTransform[n];Trtran=Transpose[Tr];
Return[Table[Trtran.Ki.Tr,{Ki,notK}]]]

GOSpBasisNoUN[N_]:=Module[{Kold},
	Kold=GOSpBasis[N];
Return[Table[If[SymmetricMatrixQ[old],old],{old,Kold}]//DeleteCases[#,Null]&]]

GOSpBasisNoU1[N_]:=Module[{Kold,notK1,notK2,Knew},
	Kold=GOSpBasis[N];
	notK1=Table[If[SymmetricMatrixQ[old]==False,old],{old,Kold}]//DeleteCases[#,Null]&;
	notK2=Table[If[Length[notK["NonzeroValues"]]==2,notK],{notK,notK1}]//DeleteCases[#,Null]&;
Return[DeleteCases[Kold,Alternatives@@notK2]]];

GOOBasis[n_]:=Module[{notK,Tr,Trtran},
	notK=Flatten[Table[SparseArray[{{i,j}->1,{j,i}->-1},{2n,2n}],{i,1,2n},{j,i+1,2n}],1];
	Tr=GOTransform[n];Trtran=Transpose[Tr];
Return[Table[Trtran.Ki.Tr,{Ki,notK}]]];

(* Matrix exp approximation *)
GOapproxExp[t_,x_]:=Module[{dim,id},
	dim=Dimensions[x][[1]];id=IdentityMatrix[dim];
Return[(id+t/2 x).Inverse[(id-t/2 x)]]];

(* Generate random Lie group elements *)
GOMranSp[n_]:=Module[{m1,m2,m3},
	m1=RandomReal[{0,1},{n,n}];
	m2=m1+Transpose[m1];
	m3=GO\[CapitalOmega]qpqp[n/2].m2;
Return[MatrixExp[m3]]];

GOMranO[n_]:=Module[{m},
	m=RandomReal[{-1,1},{n,n}];
Return[Orthogonalize[m]]];

(* Conditional logarithm *)
GOlogfunction[x_]:=If[x==0,0,Log[x]]

GOConditionalLog[x_]:=Module[{Diag,Tra,Mat},
	Diag=DiagonalMatrix[GOlogfunction/@Eigenvalues[x]];
	Tra=Eigenvectors[x]//Transpose//Chop;
	Mat=Tra.Diag.Inverse[Tra]//Chop];

(* Inverting a symplectic matrix *)
GOinvMSp[m_]:=Module[{\[CapitalOmega]}, \[CapitalOmega]=GO\[CapitalOmega]qpqp[Length[m]/2]; Return[-\[CapitalOmega].Transpose[m].\[CapitalOmega]]];


(* --------------------------------------------------------------------Optimization algorithm----------------------------------------------------------------------- *)

GOOptimize[function_, gradfunction_, J0_, M0_, {basisA_, basisB_}, {dimA_, dimB_}, tol_, steplimit_:\[Infinity]]:=

	Module[{G0, invG0, (* Initial covariance matrix *)
			KA, KB, K, (* Lie algebra basis *)
			Mold, Mnew, Eold, Enew,(* Updating function values *)
			grad, Normgrad, X, \[Epsilon], newM, GenerateM, invmetric, metric, (* Movement *)
			M0list, Elist, Normlist, diffE, diffNorm, (* Tracking values *)
			stepcount, dimM0, CorrList, order1, order2, order, keepnumber, loosenumber, donelist, (* Tracking trajectories *) 
			FinalE, FinalM, FinalElist, FinalNormlist}, (* Results *)
		
		dimM0=Length[M0]; CorrList=List[]; M0list=M0; donelist=List[]; Elist=List[];
		
		(* Generate Lie algebra basis *)
		KA=ToExpression[basisA][dimA]; KB=ToExpression[basisB][dimB];
		K=Which[
			(* Case 1: Only B is varied *)
			basisA=="None", Table[PadLeft[Ki,{Length[J0],Length[J0]}],{Ki,KB}],
			
			(* Case 2: Only A is varied *)
			basisB=="None", Table[PadRight[Ki,{Length[J0],Length[J0]}],{Ki,KA}],
			
			(* Case 3: Both A and B are varied *)
			True, Table[ArrayFlatten[{{KAi,0},{0,KBi}}],{KAi,KA},{KBi,KB}]
		];
		
		K=If[dimB==0, KA, K];
		
		(* Define new matrix generation *)
		Which[
			(* Case 1: Only B is varied *)
			basisA=="None", newM[m_,s_,x_]:=Module[{XX},
			XX=x[[2dimA+1;;Length[J0],2dimA+1;;Length[J0]]]; m.ArrayFlatten[{{IdentityMatrix[2dimA],0},{0,GOapproxExp[s,XX]}}]//SparseArray];
			
			(* Case 2: Only A is varied *)
			basisB=="None", newM[m_,s_,x_]:=Module[{XX},
			XX=x[[1;;2dimA,1;;2dimA]]; m.ArrayFlatten[{{GOapproxExp[s,XX],0},{0,IdentityMatrix[2dimB]}}]//SparseArray];
			
			(* Case 3: Both A and B are varied *)
			True, newM[m_,s_,x_]:=m.GOapproxExp[s,x]];
		
		(* Sub-routine *)
		GenerateM[\[Epsilon]_,Mold_,Mnew_,Eold_,Enew_,X_]:=Module[{s=\[Epsilon], enew=Enew, corr=0, mnew=Mnew},
			While[enew>Eold, s=s/2; corr++; mnew=newM[Mold,s,X]; enew=function[mnew,J0];]; AppendTo[CorrList,corr]; Return[{mnew//SparseArray,enew}]];
			
		(* Define natural metric *)
		G0=J0.GO\[CapitalOmega]qpqp[Length[J0]/2]; invG0=Inverse[G0];
		metric=Table[Tr[K1.K2+K1.G0.Transpose[K2].invG0+G0.Transpose[K1].invG0.K2+Transpose[K2.K1]],{K1,K},{K2,K}]; invmetric=PseudoInverse[metric];
		
		(* --------Iteration-------- *)
		
		(* Define function values and gradient for initial values *) 
		Mold=M0; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; grad=(invmetric.#)&/@grad; Normgrad=Norm/@grad; X=SparseArray/@MapThread[(-#1).K/Norm[#2]&,{grad,Normgrad}]; 
		
		Elist={Eold}//Transpose; Normlist={Normgrad}//Transpose; diffE=Eold; diffNorm=Normgrad; 
		
		(* Initialise step count *)
		stepcount=0;
	
		
		(* -----Main routine----- *)
		
		(* Stopping condition *)
		While[AllTrue[Normgrad,#>tol&] && stepcount < steplimit && AllTrue[diffE,#>10^-10&],
		
			(* Choose initial step size *)
			\[Epsilon]=Normgrad/2;
			
			(* Calculate new transformation / function value *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			{Mnew,Enew}=MapThread[GenerateM,{\[Epsilon],Mold,Mnew,Eold,Enew,X}]//Transpose;
		
		(* Define function values and gradient for new values *) 
		Mold=Mnew; Eold=Enew; grad=gradfunction[#,J0,K]&/@Mold; grad=(invmetric.#)&/@grad; Normgrad=Norm/@grad; X=SparseArray/@MapThread[(-#1).K/Norm[#2]&,{grad,Normgrad}];
		Elist=Elist//Transpose; Elist=AppendTo[Elist,Eold]//Transpose; Normlist=Normlist//Transpose; Normlist=AppendTo[Normlist,Normgrad]//Transpose;
		diffE=(Abs[#[[-1]]-#[[-2]]]/#[[-1]])&/@Elist; diffNorm=Abs[#[[-1]]-#[[-2]]]&/@Normlist;
		
		(* Update step count *)
		stepcount++;
											
		(* Checking trajectories *)
		If[stepcount/5//IntegerQ && Length[M0list]>1,
			
			(* Retain most promising 20% *)
			keepnumber=If[(.2 dimM0//Round)==0,1,.2 dimM0//Round]; 
			
			order1=Ordering[Eold,All,Less]; order2=Ordering[Normgrad,All,Greater]; 
			If[keepnumber>1,order=Join[Take[order1,keepnumber],Take[order2,keepnumber]],order=Take[order1,keepnumber]];
			
			{Mold, Eold, grad, Normgrad, X, M0list, Elist, Normlist, diffE, diffNorm}=(#[[order]])&/@{Mold, Eold, grad, Normgrad, X, M0list, Elist, Normlist, diffE, diffNorm}; dimM0=Length[order];
			];			
		];
	
	(* Check overall stopping criteria *)
	If[#<tol,AppendTo[donelist,Position[Normgrad,#]]]&/@Normgrad;
	If[stepcount >= steplimit,AppendTo[donelist,Position[Eold,Min[Eold]]];Print["Optimisation terminated by step limit."]];
	If[#<10^-10,AppendTo[donelist,Position[diffE,#]]]&/@diffE;
	
	FinalE=Min[Eold[[donelist//Flatten]]]; FinalM=M0list[[donelist//Flatten]]; FinalElist=Elist[[donelist//Flatten]]//Flatten; FinalNormlist=Normlist[[donelist//Flatten]]//Flatten; 
	
	(* --- Output --- *)

	Return[{FinalE,FinalM,CorrList,FinalElist,FinalNormlist,metric}]
];


(* --------------------------------------------------------------------Purifications and the standard form-------------------------------------------------------------------- *)

GOExtractStdForm[G_]:=Module[{NT=Dimensions[G][[1]]/2,\[CapitalOmega]T,J,RI,SWT,tra,resc,Tra,GAB,ON,TRA,check,checklist,DiagElems,Diag,Mtra,rlist}, 
	\[CapitalOmega]T=GO\[CapitalOmega]qpqp[NT]; J=-G.\[CapitalOmega]T;

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
Return[{rlist,Mtra}]];

(* Purified state for bosons *)
GOConstructPurificationBos[rlist_,{n_,m_}]:=Module[{cosh,sinh,diag,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0},
	cosh=Flatten[Table[{Cosh[2 rr],Cosh[2 rr]},{rr,rlist}]];
	sinh=Flatten[Table[{Sinh[2 rr],-Sinh[2 rr]},{rr,rlist}]];

	Q14=DiagonalMatrix[Hold/@cosh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sinh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q5=If[m==n,Null,IdentityMatrix[2(m-n)]//SparseArray];

	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{Q23,Q14}}]//SparseArray, G0=ArrayFlatten[{{Q14,Q23,0},{Q23,Q14,0},{0,0,Q5}}]//SparseArray];

	\[CapitalOmega]0=GO\[CapitalOmega]qpqp[m+n];
Return[-G0.\[CapitalOmega]0]];

(* Purified state for fermions *)
GOConstructPurificationFerm[rlist_,{n_,m_}]:=Module[{cos,sin,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0},
	cos=Table[ArrayFlatten[{{0,Cos[2 rr]},{-Cos[2 rr],0}}],{rr,rlist}];
	sin=Table[ArrayFlatten[{{0,Sin[2 rr]},{Sin[2 rr],0}}],{rr,rlist}];

	Q14=DiagonalMatrix[Hold/@cos]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sin]//ReleaseHold//ArrayFlatten//SparseArray;
	
	Q5=If[m==n,Null,ArrayFlatten[DiagonalMatrix[ConstantArray[a,m-n]]/.a -> {{0, 1},{-1, 0}}]];
	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{-Q23,Q14}}]//SparseArray, 
	G0=ArrayFlatten[{{Q14,Q23,0},{-Q23,Q14,0},{0,0,Q5}}]//SparseArray];
	
	\[CapitalOmega]0=GO\[CapitalOmega]qpqp[m+n];
Return[-G0.\[CapitalOmega]0]];


(* --------------------------------------------------------------------------Applications-------------------------------------------------------------------------- *)

(* --------------- Entanglement of Purification ---------------- *)
(* Function to generate restriction - for use in EoP *)
GORestrictionAA[partition_]:=Module[{n,n1,n2,A1,B1,A2,B2},
	{A1,B1,A2,B2}=partition;
	n=A1+B1+A2+B2; n1=A1+B1;
	Join[Table[x,{x,2A1}],Table[y+2n1,{y,2A2}]]];

(* Functionals and gradients for bosonic and fermionic EoP *)
GOEoPBos[ResAA_]:=Function[{M,J0}, Module[{n,D,logDD},
	n=Length[J0];
	D=(M.(IdentityMatrix[n]+I J0).GOinvMSp[M]/2)[[ResAA,ResAA]];
	logDD=GOConditionalLog[D.D];
	Re[1/2 Tr[D.logDD]]//Chop]];
GOEoPFerm[ResAA_]:=Function[{M,J0}, Module[{n,D,logD},
	n=Length[J0];
	D=(M.(IdentityMatrix[n]+I J0).GOinvMSp[M]/2)[[ResAA,ResAA]];
	logD=GOConditionalLog[D];
	Re[-Tr[D.logD]]//Chop]];
GOEoPgradBos[ResAA_]:=Function[{M,J0,K},Module[{n,nn,invM,dJ,D,dD,logDD},
	n=Length[J0];
	invM=GOinvMSp[M];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,K}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logDD=GOConditionalLog[D.D];
	Table[Re[1/2 Tr[dDi.logDD]]//Chop,{dDi,dD}]]];
GOEoPgradFerm[ResAA_]:=Function[{M,J0,K},Module[{n,nn,invM,KI,dJ,D,dD,logD},
	n=Length[J0]; nn=Length[J0];
	invM=GOinvMSp[M];
	KI=Table[PadLeft[Ki,{n,n}],{Ki,K}];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,KI}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logD=GOConditionalLog[D];
	Table[Re[-Tr[dDi.logD]]//Chop,{dDi,dD}]]];
	
(* --------------- Complexity of Purification ---------------- *)
(* Functionals and gradients for bosonic and fermionic CoP *)
GOCoPBos[JT_]:=Function[{M,Jref0},Module[{invM,D},
	D=M.Jref0.GOinvMSp[M].Inverse[JT]; 
	Re[Sqrt[Total[Log[#1]^2&/@Eigenvalues[D]]/8]]]];
GOCoPgradBos[JT_]:=Function[{M,Jref0,K},Module[{dimA,dimB,invM,invJT,D,invD,dD},
	dimB=Length[K[[1]]]; dimA=Length[M]-dimB;
	invM=GOinvMSp[M]; invJT=Inverse[JT];
	D=M.Jref0.invM.invJT; invD=Inverse[D];
	dD=Table[M.(KIi.Jref0-Jref0.KIi).invM.invJT,{KIi,K}];
	Table[Re[2Tr[GOConditionalLog[D].invD.dDi]],{dDi,dD}]]];

(* --------------- Energy of quadratic Hamiltonians ---------------- *)
GOenergyBos[h_]:=Function[{M,J0},Module[{\[CapitalOmega]0,G},
	\[CapitalOmega]0=GO\[CapitalOmega]qpqp[Length[M]/2];
	G=M.J0.\[CapitalOmega]0.Transpose[M]; 1/4 Tr[h.G]]];
GOenergygradBos[h_]:=Function[{M,J0,K},Module[{\[CapitalOmega]0},
	\[CapitalOmega]0=GO\[CapitalOmega]qpqp[Length[M]/2];
	Table[1/4 Tr[h.M.(KK.J0.\[CapitalOmega]0+J0.\[CapitalOmega]0.Transpose[KK]).Transpose[M]],{KK,K}]]];
GOenergyFerm[h_]:=Function[{M,J0},Module[{J},
	J=M.J0.Inverse[M]; 1/4 Tr[h.J]]];
GOenergygradFerm[h_]:=Function[{M,J0,K},
	Table[1/4 Tr[h.M.(KK.J0-J0.KK).Inverse[M]],{KK,K}]];


End[]

EndPackage[]
