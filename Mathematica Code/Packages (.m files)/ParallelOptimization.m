(* ::Package:: *)

BeginPackage["ParallelOptimization`",{"FundamentalDefinitions`"}]

ParallelOptimise::usage="Optimisation of a scalar function over the manifold of Gaussian states. 8 input arguments are (1) function to be optimised f[ M, J0], (2) associated gradient function df[ M, J0, Kbasis], (3) initial J, (4) initial M, (5) specification of basis types for subsystems A and B { basisA, basisB} (must be strings, specify 'None' for none, (6) dimensions of subsystems { dimA, dimB}, (7) tolerance on norm of gradient, and (8) step limit (default=infinity)"

Begin["Private`"]

ParallelOptimise[function_, gradfunction_, J0_, M0_, {basisA_, basisB_}, {dimA_, dimB_}, tol_, steplimit_:\[Infinity]]:=

	Module[{dim, Ktype, KA, KB, K, Mold, Mnew, Eold, Enew, FinalE, Elist, FinalElist, grad, X, \[Epsilon], stepcount, newM, 
	MnewList, NormList, FinalNormList, dimM0, ordering, StopQ, Nulls, NullsList, CorrList, DelList, pos, M0list, FinalM0list},
		
		dimM0=Length[M0]; FinalElist=List[]; FinalNormList=List[]; CorrList=List[]; M0list=M0; FinalM0list=List[];
		
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
		
		(* --------Iteration-------- *)
		
		(* Define function values and gradient for initial values *) 
		Mold=M0; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; X=(-#).K/Norm[#]&/@grad; 
		
		(* Append function value to funtion value list *)
		Elist=Partition[Eold,1]; NormList=Partition[Norm[#]&/@grad,1];
		
		(* Initialise step count *)
		stepcount=0;
		
		(* -----Main routine----- *)
		
		(* Stopping condition *)
		While[Length[Elist]>0 && stepcount < steplimit,
		(* While[Length[Elist]\[Equal]Length[M0] && stepcount<steplimit, *)
		
			(* Choose initial steo size *)
			\[Epsilon]=Norm/@grad;
			
			(* Calculate new transformation / function value *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			MnewList=List[];
			Do[
			Module[{eold=item[[1]],enew=item[[2]],mold=item[[3]],mnew=item[[4]],x=item[[5]],s=item[[6]],corr},
				corr=0;
				While[enew>eold,
					s=s/2; corr++;
					mnew=newM[mold,s,x]; enew=function[mnew,J0];
				];
				AppendTo[MnewList,mnew];
				AppendTo[CorrList,corr];
			];
			,{item,Thread[{Eold,Enew,Mold,Mnew,X,\[Epsilon]}]}
			];
		
		(* Define function values and gradient for new values *) 
		Mold=MnewList; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; X=(-#).K/Norm[#]&/@grad; 
		
		(* Append function value to funtion value list and same for norm *)
		Elist=Table[AppendTo[Elist[[i]],Eold[[i]]],{i,dimM0}]; NormList=Table[AppendTo[NormList[[i]],(Norm[#]&/@grad)[[i]]],{i,dimM0}];
		
		(* Update step count *)
		stepcount++;		
		
		(* Check individual strands for stopping criterion *)
		StopQ=If[Norm[#]>tol,0,Null]&/@grad; (* - check for stopping criterion *)		
		ordering=StopQ//Ordering; (* - rearrange *)
		Nulls=Count[StopQ,Null]; dimM0=dimM0-Nulls;

		(* - drop stopped strands *)
		Mold=Mold[[ordering]]//Drop[#,-Nulls]&; Eold=Eold[[ordering]]//Drop[#,-Nulls]&; grad=grad[[ordering]]//Drop[#,-Nulls]&; X=X[[ordering]]//Drop[#,-Nulls]&;
		M0list=M0list[[ordering]]; AppendTo[FinalM0list,Take[M0list,-Nulls]]; M0list=M0list//Drop[#,-Nulls]&;
		
		(* - add stopped strands to results and remove from Elist *)
		If[Nulls!=0, FinalElist=AppendTo[FinalElist,Take[Elist,-Nulls]]; Elist=Elist[[ordering]]//Drop[#,-Nulls]&; 
					FinalNormList=AppendTo[FinalNormList,Take[NormList,-Nulls]]; NormList=NormList[[ordering]]//Drop[#,-Nulls]&;,];
		]; 
	
	FinalElist=FinalElist//Flatten[#,1]&; FinalNormList=FinalNormList//Flatten[#,1]&;
	FinalE=Take[#,-1]&/@FinalElist//Flatten;
	
	(* Output: final function value, list of all values, number of iterations *)
	Return[{FinalE, FinalElist, Length/@FinalElist, FinalNormList, CorrList, FinalM0list//Flatten[#,1]&}]
];

End[]

EndPackage[]



