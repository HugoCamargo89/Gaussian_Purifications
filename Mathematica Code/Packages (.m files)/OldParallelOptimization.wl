(* ::Package:: *)

BeginPackage["OldParallelOptimization`",{"FundamentalDefinitions`"}]

OldParallelOptimise::usage="Optimisation of a scalar function over the manifold of Gaussian states. 8 input arguments are (1) function to be optimised f[ M, J0], (2) associated gradient function df[ M, J0, Kbasis], (3) initial J, (4) initial M, (5) specification of basis types for subsystems A and B { basisA, basisB} (must be strings, specify 'None' for none, (6) dimensions of subsystems { dimA, dimB}, (7) tolerance on norm of gradient, and (8) step limit (default=infinity)"

Begin["Private`"]

OldParallelOptimise[function_, gradfunction_, J0_, M0_, {basisA_, basisB_}, {dimA_, dimB_}, tol_, steplimit_:\[Infinity]]:=

	Module[{dim, Ktype, KA, KB, K, gradK, Mold, Mnew, Eold, Enew, Elist, grad, X, \[Epsilon], stepcount, newM, MnewList},
		
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
		
		(* Define new matrix generation *)
		Which[
			(* Case 1: Only B is varied *)
			basisA=="None", newM[m_,s_,x_]:=Module[{XX},
			XX=x[[2dimA+1;;Length[J0],2dimA+1;;Length[J0]]]; m.ArrayFlatten[{{IdentityMatrix[2dimA],0},{0,approxExp[s,XX]}}]//SparseArray];
			
			(* Case 2: Only A is varied *)
			basisB=="None", newM[m_,s_,x_]:=Module[{XX},
			XX=x[[1;;2dimA,1;;2dimA]]; m.ArrayFlatten[{{approxExp[s,XX],0},{0,IdentityMatrix[2dimB]}}]//SparseArray];
			
			(* Case 3: Both A and B are varied *)
			True, newM[m_,s_,x_]:=m.approxExp[s,x]];
		
		(* Iteration *)
		Mold=M0; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; X=(-#).K/Norm[#]&/@grad; (* Define function values and gradient for initial values *) 
		Elist={Eold}; (* Append function value to funtion value list *)
		stepcount=0; (* Initialise step count *)
		
		(* Main routine, with stopping criterion (tolerance on gradient norm / step number limit) *)
		While[!NoneTrue[Norm[#]>tol&/@grad,TrueQ] && stepcount < steplimit,

			\[Epsilon]=Norm/@grad; (* arbitrary initial choice of step size *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			MnewList=List[];
			Do[
			Module[{eold=item[[1]],enew=item[[2]],mold=item[[3]],mnew=item[[4]],x=item[[5]],s=item[[6]]},
				While[enew>eold,
					s=s/2;
					mnew=newM[mold,s,x]; enew=function[mnew,J0];
				];
				AppendTo[MnewList,mnew]
			];
			,{item,Thread[{Eold,Enew,Mold,Mnew,X,\[Epsilon]}]}
			];
		
		Mold=MnewList; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; X=(-#).K/Norm[#]&/@grad; (* Define function values and gradient for new values *) 
		AppendTo[Elist,Eold]; (* Append function value to funtion value list *)
		stepcount++ (* Update step count *)
		]; 
	
	(* Output: final function value, list of all values, number of iterations *)
	Return[{Enew, Elist, Length[Elist]}]
];

End[]

EndPackage[]



















