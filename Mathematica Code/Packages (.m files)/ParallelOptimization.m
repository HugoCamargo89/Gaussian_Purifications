(* ::Package:: *)

BeginPackage["ParallelOptimization`",{"FundamentalDefinitions`"}]

ParallelOptimise::usage="Optimisation of a scalar function over the manifold of Gaussian states. 8 input arguments are (1) function to be optimised f[ M, J0], (2) associated gradient function df[ M, J0, Kbasis], (3) initial J, (4) initial M, (5) specification of basis types for subsystems A and B { basisA, basisB} (must be strings, specify 'None' for none, (6) dimensions of subsystems { dimA, dimB}, (7) tolerance on norm of gradient, and (8) step limit (default=infinity)"

Begin["Private`"]

ParallelOptimise[function_, gradfunction_, J0_, M0_, {basisA_, basisB_}, {dimA_, dimB_}, tol_, steplimit_:\[Infinity]]:=

	Module[{KA, KB, K, (* Lie algebra basis *)
			Mold, Mnew, Eold, Enew,(* Updating function values *)
			grad, Normgrad, X, \[Epsilon], newM, GenerateM, (* Movement *)
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
		
		Elist={Eold}//Transpose; Normlist={Normgrad}//Transpose; diffE=Eold; diffNorm=Normgrad; 
		
		(* Initialise step count *)
		stepcount=0;
	
		
		(* -----Main routine----- *)
		
		(* Stopping condition *)
		While[AllTrue[Normgrad,#>tol&] && stepcount < steplimit && AllTrue[diffE,#>10^-10&],
		
			(* Choose initial step size *)
			\[Epsilon]=Normgrad;
			
			(* Calculate new transformation / function value *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			{Mnew,Enew}=MapThread[GenerateM,{\[Epsilon],Mold,Mnew,Eold,Enew,X}]//Transpose;
		
		(* Define function values and gradient for new values *) 
		Mold=Mnew; Eold=Enew; grad=gradfunction[#,J0,K]&/@Mold; Normgrad=Norm/@grad; X=SparseArray/@MapThread[(-#1).K/Norm[#2]&,{grad,Normgrad}];
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
	
	(* Check overall stopping criterion *)
	If[#<tol,AppendTo[donelist,Position[Normgrad,#]]]&/@Normgrad;
	If[stepcount >= steplimit,AppendTo[donelist,Position[Eold,Min[Eold]]];Print["Optimisation terminated by step limit."]];
	If[#<10^-10,AppendTo[donelist,Position[diffE,#]]]&/@diffE;
	
	FinalE=Min[Eold[[donelist//Flatten]]]; FinalM=M0list[[donelist//Flatten]]; FinalElist=Elist[[donelist//Flatten]]//Flatten; FinalNormlist=Normlist[[donelist//Flatten]]//Flatten; 
	
	(* --- Output --- *)

	Return[{FinalE,FinalM,CorrList//Total,FinalElist,FinalNormlist}]
];

End[]

EndPackage[]



