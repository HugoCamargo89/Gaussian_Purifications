(* ::Package:: *)

BeginPackage["ParallelOptimization`",{"FundamentalDefinitions`"}]

ParallelOptimise::usage="Optimisation of a scalar function over the manifold of Gaussian states. 8 input arguments are (1) function to be optimised f[ M, J0], (2) associated gradient function df[ M, J0, Kbasis], (3) initial J, (4) initial M, (5) specification of basis types for subsystems A and B { basisA, basisB} (must be strings, specify 'None' for none, (6) dimensions of subsystems { dimA, dimB}, (7) tolerance on norm of gradient, and (8) step limit (default=infinity)"

Begin["Private`"]

ParallelOptimise[function_, gradfunction_, J0_, M0_, {basisA_, basisB_}, {dimA_, dimB_}, tol_, steplimit_:\[Infinity]]:=

	Module[{dim, Ktype, KA, KB, K, (* Lie algebra basis *)
			Mold, Mnew, Eold, Enew, (* Updating function values *)
			grad, Normgrad, X, \[Epsilon], newM, GenerateM, (* Movement *)
			Elist, FinalElist, MnewList, NormList, FinalNormList, M0list, FinalM0list, (* Tracking values *)
			stepcount, dimM0, initialdimM0, CorrList, order, keepnumber, loosenumber}, (* Tracking trajectories *) 
		
		initialdimM0=Length[M0]; dimM0=Length[M0]; CorrList=List[]; M0list=M0; FinalElist=List[]; FinalNormList=List[];
		
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
		
		(* Sub-routine *)
		GenerateM[\[Epsilon]_,Mold_,Mnew_,Eold_,Enew_,X_]:=Module[{s=\[Epsilon], enew=Enew, corr=0, mnew=Mnew},
			While[enew>Eold, s=s/2; corr++; mnew=newM[Mold,s,X]; enew=function[mnew,J0];]; AppendTo[CorrList,corr]; Return[{mnew//SparseArray,enew}]];
		
		
		(* --------Iteration-------- *)
		
		(* Define function values and gradient for initial values *) 
		Mold=M0; Eold=function[#,J0]&/@Mold; grad=gradfunction[#,J0,K]&/@Mold; Normgrad=Norm/@grad; X=SparseArray/@MapThread[(-#1).K/Norm[#2]&,{grad,Normgrad}]; 
		
		(* Append function value to funtion value list *)
		Elist=Partition[Eold,1]; NormList=Partition[Normgrad,1];
		
		(* Initialise step count *)
		stepcount=0;
		
		(* -----Main routine----- *)
		
		(* Stopping condition *)
		While[AllTrue[Normgrad,#>tol&] && stepcount < steplimit,
		
			(* Choose initial step size *)
			\[Epsilon]=Normgrad;
			
			(* Calculate new transformation / function value *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			{Mnew,Enew}=MapThread[GenerateM,{\[Epsilon],Mold,Mnew,Eold,Enew,X}]//Transpose;
		
		(* Define function values and gradient for new values *) 
		Mold=Mnew; Eold=Enew; grad=gradfunction[#,J0,K]&/@Mold; Normgrad=Norm/@grad; X=SparseArray/@MapThread[(-#1).K/Norm[#2]&,{grad,Normgrad}];
		
		(* Append function value to funtion value list and same for norm *)
		Elist=Table[AppendTo[Elist[[i]],Eold[[i]]],{i,dimM0}]; NormList=Table[AppendTo[NormList[[i]],Normgrad[[i]]],{i,dimM0}];
		
		(* Update step count *)
		stepcount++;		
		
		(* Checking trajectories *)
		If[stepcount/4//IntegerQ && Length[Elist]>1,
			
			(* Retain most promising 20% *)
			keepnumber=If[(.2 dimM0//Round)==0,1,.2 dimM0//Round]; loosenumber=dimM0-keepnumber;
			
			order=Ordering[Eold,All,Less];
			FinalElist=AppendTo[FinalElist,Take[Elist[[order]],-loosenumber]]; FinalNormList=AppendTo[FinalNormList,Take[NormList[[order]],-loosenumber]];
			{Mold, Eold, grad, Normgrad, X, M0list}=Drop[#[[order]],-loosenumber]&/@{Mold, Eold, grad, Normgrad, X, M0list}; dimM0=dimM0-loosenumber;
			];	
		];
		If[Length[Elist]==1, FinalElist=AppendTo[FinalElist,Take[Elist,1]]; FinalNormList=AppendTo[FinalNormList,Take[NormList,1]]; ];
	
	FinalElist=FinalElist//Flatten[#,1]&; FinalNormList=FinalNormList//Flatten[#,1]&;
	
	(* --- Output --- *)
	(* Return[{Take[#,-1]&/@FinalElist//Flatten, FinalElist, Length/@FinalElist, CorrList, FinalNormList, Flatten[FinalM0list,1]}] *)
	Return[{Eold[[All]], FinalElist, Length/@FinalElist, FinalNormList, CorrList, M0list[[1]]}]
];

End[]

EndPackage[]



