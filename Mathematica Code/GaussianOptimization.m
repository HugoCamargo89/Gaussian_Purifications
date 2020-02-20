(* ::Package:: *)

BeginPackage["GaussianOptimization`"];


(* ::Text:: *)
(*Usage definitions*)


(* Optimization algorithm *)
GOOptimize::usage="The function GOOptimize[ProblemSpecific,SystemSpecific,ProcedureSpecific] provides a general gradient descent optimization algorithm and sits at the core of the GO package.

{FinalF,FinalM,CorrList,FinalFlist,Flists,FinalNormlist}=GOOptimize[ProblemSpecific,SystemSpecific,ProcedureSpecific]

\!\(\*
StyleBox[\"The\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"input\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"consists\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"of\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"the\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"following\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"pieces\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\":\",\nFontWeight->\"Bold\"]\)
ProblemSpecific={\!\(\*
StyleBox[\"F\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"dF\",\nFontWeight->\"Bold\"]\)}
	\!\(\*
StyleBox[\"F\",\nFontWeight->\"Bold\"]\) represents the function to be minimized
	\!\(\*
StyleBox[\"dF\",\nFontWeight->\"Bold\"]\) represents the gradient function of f
SystemSpecific={\!\(\*
StyleBox[\"J0\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"listM0\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"KK\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"geometry\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"newM\",\nFontWeight->\"Bold\"]\)}
	\!\(\*
StyleBox[\"J0\",\nFontWeight->\"Bold\"]\) represents a reference point
	\!\(\*
StyleBox[\"listM0\",\nFontWeight->\"Bold\"]\) represents a list of starting points, i.e., of transformations from the reference point J0
	\!\(\*
StyleBox[\"KK\",\nFontWeight->\"Bold\"]\) represents a local basis of tangent space (e.g. subspace of Lie algebra)
	\!\(\*
StyleBox[\"geometry\",\nFontWeight->\"Bold\"]\) represents a local notion of geometry (e.g. metric) with respect to KK
	\!\(\*
StyleBox[\"newM\",\nFontWeight->\"Bold\"]\) represents the function that takes the gradient as input and constructs a new point on the state manifold
ProcedureSpecific={\!\(\*
StyleBox[\"tolF\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"tolgrad\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"steplimit\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"stepsizefunction\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"trackall\",\nFontWeight->\"Bold\"]\)}
	\!\(\*
StyleBox[\"tolF\",\nFontWeight->\"Bold\"]\) represents the stopping tolerance of the function value of F
	\!\(\*
StyleBox[\"tolgrad\",\nFontWeight->\"Bold\"]\) represents the stopping tolerance of the gradient norm of dF
	\!\(\*
StyleBox[\"steplimit\",\nFontWeight->\"Bold\"]\) represents a maximal number of steps (typically, we choose \[Infinity])
	\!\(\*
StyleBox[\"trackall\",\nFontWeight->\"Bold\"]\) is an option to track all trajectories starting from listM0 instead of dropping some

\!\(\*
StyleBox[\"The\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"output\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"consists\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"of\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"the\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"following\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"pieces\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\":\",\nFontWeight->\"Bold\"]\)
{\!\(\*
StyleBox[\"FinalF\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"FinalM\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"CorrList\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"FinalFlist\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"Flists\",\nFontWeight->\"Bold\"]\),\!\(\*
StyleBox[\"FinalNormlist\",\nFontWeight->\"Bold\"]\)}
	\!\(\*
StyleBox[\"FinalF\",\nFontWeight->\"Bold\"]\) represents the final function value of F
	\!\(\*
StyleBox[\"FinalM\",\nFontWeight->\"Bold\"]\) represents the final transformation from J0 to the final state
	\!\(\*
StyleBox[\"CorrList\",\nFontWeight->\"Bold\"]\) lists the number of stepsize corrections at each step of the final trajectory
	\!\(\*
StyleBox[\"FinalFlists\",\nFontWeight->\"Bold\"]\) lists the successive function values in the optimal trajectory
	\!\(\*
StyleBox[\"Flists\",\nFontWeight->\"Bold\"]\) lists all successive function values for all trajectories (Lucas: is that right?)
	\!\(\*
StyleBox[\"FinalNormlist\",\nFontWeight->\"Bold\"]\) lists the successive gradient norms in the optimal trajectory
";


(* Symplectic forms / Covariance matrices *)
GO\[CapitalOmega]qqpp::usage="GO\[CapitalOmega]qqpp[N] generates \[CapitalOmega] in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GO\[CapitalOmega]qpqp::usage="GO\[CapitalOmega]qpqp[N] generates \[CapitalOmega] in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GO\[CapitalOmega]aabb::usage="GO\[CapitalOmega]aabb[N] generates \[CapitalOmega] in basis (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GO\[CapitalOmega]abab::usage="GO\[CapitalOmega]abab[N] generates \[CapitalOmega] in basis (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";

GOGqqpp::usage="GOGqqpp[N] generates G in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOGqpqp::usage="GOGqpqp[N] generates G in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOGaabb::usage="GOGaabb[N] generates G in basis (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GOGabab::usage="GOGabab[N] generates G in basis (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";


(* Basis Transformation matrices *)
(* Re-ordering *)
GOqqppFROMqpqp::usage="GOqqppFROMqpqp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOqpqpFROMqqpp::usage="GOqpqpFROMqqpp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOaabbFROMabab::usage="GOaabbFROMabab[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GOababFROMaabb::usage="GOababFROMaabb[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
(* q to a *)
GOababFROMqpqp::usage="GOababFROMqpqp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GOaabbFROMqpqp::usage="GOaabbFROMqpqp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GOababFROMqqpp::usage="GOababFROMqqpp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
GOaabbFROMqqpp::usage="GOaabbFROMqqpp[N] generates a basis transformation (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...) for N deg. of freedom";
(* a to q *)
GOqpqpFROMabab::usage="GOqpqpFROMabab[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOqqppFROMabab::usage="GOqqppFROMabab[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOqpqpFROMaabb::usage="GOqpqpFROMaabb[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";
GOqqppFROMaabb::usage="GOqqppFROMaabb[N] generates a basis transformation (\!\(\*SubscriptBox[\(a\), \(1\)]\),\!\(\*SubscriptBox[\(a\), \(2\)]\),...,\!\(\*SubscriptBox[\(b\), \(1\)]\),\!\(\*SubscriptBox[\(b\), \(2\)]\),...)\[Rule](\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),...,\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...) for N deg. of freedom";


(* Transforming between structures *)
GOTransformGtoJ::usage="GOTransformGtoJ[G,'Gbasis','Jbasis'] generates J in Jbasis associated with G in Gbasis";
GOTransform\[CapitalOmega]toJ::usage="GOTransform\[CapitalOmega]toJ[\[CapitalOmega],'\[CapitalOmega]basis','Jbasis'] generates J in Jbasis associated with \[CapitalOmega] in \[CapitalOmega]basis";
GOTransformJtoG::usage="GOTransformJtoG[J,'Jbasis','Gbasis'] generates G in Gbasis associated with J in Jbasis";
GOTransformJto\[CapitalOmega]::usage="GOTransformJto\[CapitalOmega][J,'Jbasis','\[CapitalOmega]basis'] generates \[CapitalOmega] in \[CapitalOmega]basis associated with J in Jbasis";


(* Random transformations *)
GORandomTransformation::usage="GORandomTransformation[n,'Basis'] generates a random transformation in the Lie group with algebra 'Basis'";


(* Tools for computational efficiency *)
GOapproxExp::usage="GOapproxExp[\[Epsilon],K]=(1 + \[Epsilon]/2 K)/(1 - \[Epsilon]/2 K) approximates MatrixExponential[\[Epsilon]X]";
GOlogfunction::usage="GOlogfunction[x] takes value 0 when x=0 and log[x] otherwise";
GOConditionalLog::usage="GOConditionalLog[x] calculates MatrixLog[x] by applying logfunction to the eigenvalues of x";
GOinvMSp::usage="GOinvMSp[m,'Mform'] inverts a symplectic matrix m in the basis Mform as \!\(\*SuperscriptBox[\(m\), \(-1\)]\)=-\[CapitalOmega]m\[CapitalOmega]";


(* Lie algebra bases *)
GOLieBasisSp::usage="GOLieBasisSp[N] generates Sp(2N,R) Lie algebra basis in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)";
GOLieBasisSpNoUN::usage="GOLieBasisSpNoUN[N] generates Sp(2N,R)/U(N) Lie algebra basis in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)";
GOLieBasisSpNoU1::usage="GOLieBasisSpNoU1[N] generates Sp(2N,R)/U(1\!\(\*SuperscriptBox[\()\), \(N\)]\) Lie algebra basis in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)";
GOLieBasisSpMixing::usage="GOLieBasisSpMixing[\!\(\*SubscriptBox[\(N\), \(A\)]\),\!\(\*SubscriptBox[\(N\), \(B\)]\)] generates the compound Lie algebra basis sp(2(\!\(\*SubscriptBox[\(N\), \(A\)]\)+\!\(\*SubscriptBox[\(N\), \(B\)]\)),R)/sp(2\!\(\*SubscriptBox[\(N\), \(A\)]\),R) x sp(2\!\(\*SubscriptBox[\(N\), \(B\)]\),R)";
GOLieBasisSO::usage="GOLieBasisSO[N] generates O(2N,R) Lie algebra basis in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)";
GOLieBasisSONoUN::usage="GOLieBasisSpNoUN[N] generates SO(2N,R)/U(N) Lie algebra basis in basis (\!\(\*SubscriptBox[\(q\), \(1\)]\),\!\(\*SubscriptBox[\(p\), \(1\)]\),\!\(\*SubscriptBox[\(q\), \(2\)]\),\!\(\*SubscriptBox[\(p\), \(2\)]\),...)";
GOLieBasisEmpty::usage="GOLieBasisEmpty[N] constructs the empty Lie basis given by a list containing a single 2N by 2N matrix with zero entries. It can be used to construct a compound basis using GOLieBasisCompound.";


GOLieBasisCompound::usage="GOCompoundBasis[{LieAlgebra1,LieAlgebra2,...,LieAlgebraN}] creates the space spanned by different Lie algebra spaces collected in a vector {LieAlgebra1,LieAlgebra2,...,LieAlgebraN}. For degrees of freedom without changes, we should include GOLieAlgebraEmpty[n].";


GOMetricSp::usage="GOMetricSp[K,\!\(\*SubscriptBox[\(J\), \(0\)]\)] generates the natural metric for the symplectic manifold with Lie basis K, for an initial state \!\(\*SubscriptBox[\(J\), \(0\)]\)";
GOMetricO::usage="GOMetricO[K,\!\(\*SubscriptBox[\(J\), \(0\)]\)] generates the natural metric for the orthogonal manifold with Lie basis K";
GOGeometryConst::usage="GeometryConst[LieBasis,J0,pm] creates the list geometry={metric,invmetric} of a constant geometry around J0. It works for both, bosonic and fermionic systems by choosing pm=+1 for bosons and pm=-1 for fermions.";


(* Purifications and the standard form *)
GOExtractStdFormG::usage="{rlist,Tra}=GOExtractStdFormG[G]
Returns the parameters for constructing the standard form of G (list of the \!\(\*SubscriptBox[\(r\), \(i\)]\) and the transformation matrix Tra that puts G into standard form";
GOExtractStdForm\[CapitalOmega]::usage="{rlist,Tra}=GOExtractStdForm\[CapitalOmega][\[CapitalOmega]]
Returns the parameters for constructing the standard form of \[CapitalOmega] (list of the \!\(\*SubscriptBox[\(r\), \(i\)]\) and the transformation matrix Tra that puts \[CapitalOmega] into standard form";
GOPurifyStandardGBoson::usage="\!\(\*SubscriptBox[\(G\), \(sta\)]\)=GOPurifyStandardGBoson[rlist,N,'Gform']
Constructs the bosonic purified state convariance matrix in the basis Gform with N degrees of freedom in the ancillary";
GOPurifyStandardJBoson::usage="\!\(\*SubscriptBox[\(J\), \(sta\)]\)=GOPurifyStandardJBoson[rlist,N,'Jform']
Constructs the bosonic purified state complex structure in the basis Jform with N degrees of freedom in the ancillary";
GOPurifyStandard\[CapitalOmega]Fermion::usage="\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(sta\)]\)=GOPurifyStandard\[CapitalOmega]Fermion[rlist,N,'\[CapitalOmega]form']
Constructs the fermionic purified state convariance matrix in the basis \[CapitalOmega]form with N degrees of freedom in the ancillary";
GOPurifyStandardJFermion::usage="\!\(\*SubscriptBox[\(J\), \(sta\)]\)=GOPurifyStandardJFermion[rlist,N,'Jform']
Constructs the fermionic purified state complex structure in the basis Jform with N degrees of freedom in the ancillary";


(* Application: Entanglement of purification *)
GORestrictionAA::usage="GORestrictionAA[{dimA1, dimB1, dimA2, dimB2}] generates the restriction function for use in calculating EoP for a system decomposed into degrees of freedom {dimA1, dimB1, dimA2, dimB2}";
GOEoPBos::usage="GOEoPBos[Restriction] generates a function f[M,\!\(\*SubscriptBox[\(J\), \(0\)]\)] which calculates the bosonic EoP for the given restriction";
GOEoPFerm::usage="GOEoPBos[Restriction] generates a function f[M,\!\(\*SubscriptBox[\(J\), \(0\)]\)] which calculates the fermionic EoP for the given restriction";
GOEoPgradBos::usage="GOEoPgradBos[Restriction] generates a function df[M,\!\(\*SubscriptBox[\(J\), \(0\)]\),K,\!\(\*SuperscriptBox[\(G\), \(-1\)]\)] which calculates the gradient of the bosonic EoP for the given restriction with respect to the Lie algebra basis K";
GOEoPgradFerm::usage="GOEoPgradFerm[Restriction] generates a function df[M,\!\(\*SubscriptBox[\(J\), \(0\)]\),K,\!\(\*SuperscriptBox[\(G\), \(-1\)]\)] which calculates the gradient of the fermionic EoP for the given restriction with respect to the Lie algebra basis K";
(* Application: Complexity of purification *)
GOCoPBos::usage="GOCoPBos[\!\(\*SubscriptBox[\(J\),\(T\)]\)] generates a function f[M,\!\(\*SubscriptBox[\(J\), \(R0\)]\)] which calculates the bosonic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\)";
GOCoPgradBos::usage="GOCoPBos[\!\(\*SubscriptBox[\(J\), \(T\)]\)] generates a function df[M,\!\(\*SubscriptBox[\(J\), \(R0\)]\),K,\!\(\*SuperscriptBox[\(G\), \(-1\)]\)] which calculates the gradient of the bosonic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\) and the Lie algebra basis K";
GOCoPFerm::usage="GOCoPFerm[\!\(\*SubscriptBox[\(J\),\(T\)]\)] generates a function f[M,\!\(\*SubscriptBox[\(J\), \(R0\)]\)] which calculates the fermionic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\)";
GOCoPgradFerm::usage="GOCoPgradFerm[\!\(\*SubscriptBox[\(J\), \(T\)]\)] generates a function df[M,\!\(\*SubscriptBox[\(J\), \(R0\)]\),K,\!\(\*SuperscriptBox[\(G\), \(-1\)]\)] which calculates the gradient of the fermionic CoP with respect to the given target state \!\(\*SubscriptBox[\(J\), \(T\)]\) and the Lie algebra basis K";
(* Appliation: Energy of quadratic Hamiltonians *)
GOenergyBos::usage="GOenergyBos[h] generates an energy function f[M,\!\(\*SubscriptBox[\(J\), \(0\)]\)] for a bosonic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\)";
GOenergygradBos::usage="GOenergygradBos[h] generates an energy gradient function df[M,\!\(\*SubscriptBox[\(J\), \(0\)]\),K] for a bosonic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\) with respect to the Lie algebra basis K";
GOenergyFerm::usage="GOenergyFerm[h] generates an energy function f[M,\!\(\*SubscriptBox[\(J\), \(0\)]\)] for a fermionic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\)";
GOenergygradFerm::usage="GOenergygradBos[h] generates an energy gradient function df[M,\!\(\*SubscriptBox[\(J\), \(0\)]\),K] for a fermionic quadratic Hamiltonian H=\!\(\*SubscriptBox[\(h\), \(ab\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(a\)]\)\!\(\*SuperscriptBox[\(\[Xi]\), \(b\)]\) with respect to the Lie algebra basis K";


(* ::Text:: *)
(*Functions*)


Begin["Private`"];


(* Symplectic forms / Covariance matrices *)
GO\[CapitalOmega]qqpp[n_]:=ArrayFlatten[{{0,IdentityMatrix[n]},{-IdentityMatrix[n],0}}]//SparseArray; 
GO\[CapitalOmega]qpqp[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> ({{0, 1},{-1, 0}})]//SparseArray ;
GO\[CapitalOmega]aabb[n_]:=ArrayFlatten[{{0,-I IdentityMatrix[n]},{I IdentityMatrix[n],0}}]//SparseArray; 
GO\[CapitalOmega]abab[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> ({{0, -I},{I, 0}})]//SparseArray ;

GOGqqpp[n_]:=IdentityMatrix[2n]//SparseArray;
GOGqpqp[n_]:=IdentityMatrix[2n]//SparseArray;
GOGaabb[n_]:=ArrayFlatten[{{0,IdentityMatrix[n]},{IdentityMatrix[n],0}}]//SparseArray ;
GOGabab[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> ({{0, 1},{1, 0}})]//SparseArray ;


(* Basis Transformation matrices *)
(* Re-ordering *)
GOqqppFROMqpqp[n_]:=SparseArray[Join[Table[{k,2k-1}->1,{k,n}],Table[{n+k,2k}->1,{k,n}]]];
GOqpqpFROMqqpp[n_]:=SparseArray[Join[Table[{2k-1,k}->1,{k,n}],Table[{2k,n+k}->1,{k,n}]]];
GOaabbFROMabab[n_]:=SparseArray[Join[Table[{k,2k-1}->1,{k,n}],Table[{n+k,2k}->1,{k,n}]]];
GOababFROMaabb[n_]:=SparseArray[Join[Table[{2k-1,k}->1,{k,n}],Table[{2k,n+k}->1,{k,n}]]];
(* q to a *)
GOababFROMqpqp[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> (1/Sqrt[2]{{1, I},{1, -I}})]//SparseArray ;
GOaabbFROMqpqp[n_]:=GOqqppFROMqpqp[n].GOababFROMqpqp[n]//SparseArray ;
GOababFROMqqpp[n_]:=GOababFROMqpqp[n].GOqpqpFROMqqpp[n]//SparseArray ;
GOaabbFROMqqpp[n_]:=GOaabbFROMabab[n].GOababFROMqqpp[n]//SparseArray ;
(* a to q *)
GOqpqpFROMabab[n_]:=ArrayFlatten[DiagonalMatrix[ConstantArray[a,n]]/.a -> (1/Sqrt[2]{{1, 1},{-I, I}})]//SparseArray ;
GOqqppFROMabab[n_]:=GOqqppFROMqpqp[n].GOqpqpFROMabab[n]//SparseArray ;
GOqpqpFROMaabb[n_]:=GOqpqpFROMabab[n].GOababFROMaabb[n]//SparseArray ;
GOqqppFROMaabb[n_]:=GOqqppFROMabab[n].GOababFROMaabb[n]//SparseArray ;


(* J from G and \[CapitalOmega] *)
GOTransformGtoJ[G_,Gform_,Jform_]:=Module[{n,Mtra,Gcorrect,\[CapitalOmega]correct},
	n=Length[G]/2;
	Mtra=If[Gform==Jform,IdentityMatrix[2n],ToExpression["GO"<>Jform<>"FROM"<>Gform][n]];
	Gcorrect=Mtra.G.Transpose[Mtra];
	\[CapitalOmega]correct=ToExpression["GO\[CapitalOmega]"<>Jform][n];
	-Gcorrect.\[CapitalOmega]correct];
GOTransform\[CapitalOmega]toJ[\[CapitalOmega]_,\[CapitalOmega]form_,Jform_]:=Module[{n,Mtra,Gcorrect,\[CapitalOmega]correct},
	n=Length[\[CapitalOmega]]/2;
	Mtra=If[\[CapitalOmega]form==Jform,IdentityMatrix[2n],ToExpression["GO"<>Jform<>"FROM"<>\[CapitalOmega]form][n]];
	\[CapitalOmega]correct=Mtra.\[CapitalOmega].Transpose[Mtra];
	Gcorrect=ToExpression["GOG"<>Jform][n];
	-Gcorrect.\[CapitalOmega]correct];
GOTransformJtoG[J_,Jform_,Gform_]:=Module[{n,Mtra,Jcorrect,\[CapitalOmega]correct},
	n=Length[J]/2;
	Mtra=If[Gform==Jform,IdentityMatrix[2n],ToExpression["GO"<>Gform<>"FROM"<>Jform][n]];
	Jcorrect=Mtra.J.Transpose[Mtra];
	\[CapitalOmega]correct=ToExpression["GO\[CapitalOmega]"<>Gform][n];
	Jcorrect.\[CapitalOmega]correct];
GOTransformJto\[CapitalOmega][J_,Jform_,\[CapitalOmega]form_]:=Module[{n,Mtra,Jcorrect,Gcorrect},
	n=Length[J]/2;
	Mtra=If[\[CapitalOmega]form==Jform,IdentityMatrix[2n],ToExpression["GO"<>\[CapitalOmega]form<>"FROM"<>Jform][n]];
	Jcorrect=Mtra.J.Transpose[Mtra];
	Gcorrect=ToExpression["GOG"<>\[CapitalOmega]form][n];
	Jcorrect.Gcorrect];


(* Random transformations *)
GORandomTransformation[{dim1_,dim2_},LieBasis_]:=Module[{K,coeffs,KK,M},
	K=ToExpression[LieBasis][dim1,dim2];
	coeffs=RandomReal[{-1,1},Length[K]];
	KK=coeffs.K;
	M=MatrixExp[KK]];


(* Matrix exp approximation *)
GOapproxExp[t_,x_]:=Module[{dim,id},
	dim=Dimensions[x][[1]];id=IdentityMatrix[dim];
	Return[(id+t/2 x).Inverse[(id-t/2 x)]]];
	
(* Conditional logarithm *)
GOlogfunction[x_]:=If[x==0,0,Log[x]]
GOConditionalLog[x_]:=Module[{Diag,Tra,Mat},
	Diag=DiagonalMatrix[GOlogfunction/@Eigenvalues[x]];
	Tra=Eigenvectors[x]//Transpose//Chop;
	Mat=Tra.Diag.Inverse[Tra]//Chop];

(* Inverting a symplectic matrix *)
GOinvMSp[m_,Mform_:"qpqp"]:=Module[{\[CapitalOmega]}, \[CapitalOmega]=ToExpression["GO\[CapitalOmega]"<>Mform][Length[m]/2]; Return[-\[CapitalOmega].Transpose[m].\[CapitalOmega]]];


(* Lie algebra bases *)

(* Generates the Lie algebra of Sp(2n,R) *)
GOLieBasisSp[n_]:=Module[{SYM,ASYM,notK,Tr,Trtran},
	SYM=Flatten[Table[SparseArray[{{i,j}->If[i!=j,N[1/Sqrt[2]],1],{j,i}->If[i!=j,N[1/Sqrt[2]],1]},{n,n}],{i,1,n},{j,i,n}],1];
	ASYM=Flatten[Table[SparseArray[{{i,j}->N[1/Sqrt[2]],{j,i}->-N[1/Sqrt[2]]},{n,n}],{i,1,n},{j,i+1,n}],1];
	notK=Join[
		Table[ArrayFlatten[{{m,0},{0,-m}}],{m,SYM}],
		Table[ArrayFlatten[{{m,0},{0,m}}],{m,ASYM}],
		Table[ArrayFlatten[{{0,m},{m,0}}],{m,SYM}],
		Table[ArrayFlatten[{{0,m},{-m,0}}],{m,SYM}]];
	Tr=GOqpqpFROMqqpp[n];
Return[Table[Tr.Ki.Transpose[Tr],{Ki,notK}]]]

(* Generates the subspace of the Lie algebra of Sp(2n,R) where we remove the u(n) sub algebra*)
GOLieBasisSpNoUN[n_]:=Module[{Kold},
	Kold=GOLieBasisSp[n];
Return[Table[If[SymmetricMatrixQ[old],old],{old,Kold}]//DeleteCases[#,Null]&]];

(* Generates the subspace of the Lie algebra of Sp(2n,R) where we remove u(1)^n *)
GOLieBasisSpNoU1[n_]:=Module[{N,Kold,notK1,notK2,Knew},
	Kold=GOLieBasisSp[n];
	notK1=Table[If[SymmetricMatrixQ[old]==False,old],{old,Kold}]//DeleteCases[#,Null]&;
	notK2=Table[If[Length[notK["NonzeroValues"]]==2,notK],{notK,notK1}]//DeleteCases[#,Null]&;
Return[DeleteCases[Kold,Alternatives@@notK2]]];

(* Generates the subspace of the Lie algebra of Sp(2dim1+2dim2,R) where we remove the subalgebras Sp(2dim1,R) and Sp(2dim2,R) *)
GOLieBasisSpMixing[dim1_,dim2_]:=Module[{FullBasis,Indices1,Indices2,Indices,MQ2,MQ3,CheckM,CheckList},
	FullBasis=GOLieBasisSp[dim1+dim2];
	Indices1=Range[2dim1];Indices2=Range[2dim1+1,2(dim1+dim2)];

	MQ2[M_]:=M[[Indices1,Indices2]];MQ3[M_]:=M[[Indices2,Indices1]];
	CheckM[M_]:=Module[{Q2,Q3},
		Q2=MQ2[M]; Q3=MQ3[M]; AllTrue[Q2//Flatten,#==0&]&&AllTrue[Q2//Flatten,#==0&]];

	CheckList=CheckM/@FullBasis; Indices=Position[CheckList,False]//Flatten;
	FullBasis[[Indices]]];

(* Generates the Lie algebra of SO(2n) *)
GOLieBasisSO[n_]:=Module[{notK,Tr,Trtran},
	notK=Flatten[Table[SparseArray[{{i,j}->1,{j,i}->-1},{2n,2n}],{i,1,2n},{j,i+1,2n}],1];
	Tr=GOqpqpFROMqqpp[n];
Return[Table[Tr.Ki.Transpose[Tr],{Ki,notK}]]];

(* Generates the subspace of the Lie algebra of SO(2n) where we remove the u(n) sub algebra*)
GOLieBasisSONoUN[N_]:=Module[{FullBasis,\[CapitalOmega]},
	FullBasis=GOLieBasisSO[N];
	\[CapitalOmega]=GO\[CapitalOmega]qpqp[N];
Return[Complement[FullBasis,Select[FullBasis,SymmetricMatrixQ[#.\[CapitalOmega]]&]]]];

(* Generates the empty Lie algebra of n degrees of freedom, i.e, a list containing a single zero 2n by 2n matrix *)
GOLieBasisEmpty[n_]:={ConstantArray[0,{2n,2n}]};

(* Generating compound basis *)
GOLieBasisCompound[VecOfLieAlgebra_]:=Module[{NumberOfLieAlgebras,VecOfDim,VecOfEnlarged,TotalDim,Ktotal,FirstDim},
	NumberOfLieAlgebras=Length[VecOfLieAlgebra];
	VecOfDim=Length[#[[1]]]&/@VecOfLieAlgebra;
	TotalDim=Total[VecOfDim];
	
	VecOfEnlarged=Table[
	Which[
	
	(* Case 1: Only one component (kind of pointless, but ok)*)
	NumberOfLieAlgebras==1,VecOfLieAlgebra[1],
	
	(* Case 2: First, Lie algebra*)
	i==1, If[VecOfLieAlgebra[[i,1]]==ConstantArray[0,{VecOfDim[[i]],VecOfDim[[i]]}],{},
			Table[ArrayFlatten[({
 {Ki, 0},
 {0, ConstantArray[0,{TotalDim-VecOfDim[[i]],TotalDim-VecOfDim[[i]]}]}
})],{Ki,VecOfLieAlgebra[[i]]}]
			],
	i==NumberOfLieAlgebras,If[VecOfLieAlgebra[[i,1]]==ConstantArray[0,{VecOfDim[[i]],VecOfDim[[i]]}],{},
			Table[ArrayFlatten[({
 {ConstantArray[0,{TotalDim-VecOfDim[[i]],TotalDim-VecOfDim[[i]]}], 0},
 {0, Ki}
})],{Ki,VecOfLieAlgebra[[i]]}]
			],
	True,If[VecOfLieAlgebra[[i,1]]==ConstantArray[0,{VecOfDim[[i]],VecOfDim[[i]]}],{},
			FirstDim=Total[VecOfDim[[1;;i-1]]];
			Table[ArrayFlatten[({
 {ConstantArray[0,{FirstDim,FirstDim}], 0, 0},
 {0, Ki, 0},
 {0, 0, ConstantArray[0,{TotalDim-VecOfDim[[i]]-FirstDim,TotalDim-VecOfDim[[i]]-FirstDim}]}
})],{Ki,VecOfLieAlgebra[[i]]}]
			]
		
	],{i,1,NumberOfLieAlgebras}];
	

	Ktotal=SparseArray/@Flatten[VecOfEnlarged,1];
	Return[Ktotal]];


(* Natural metrics *)
GOMetricSp[LieBasis_,J0_]:=Module[{G0,invG0},
	G0=J0.GO\[CapitalOmega]qpqp[Length[J0]/2]; invG0=Inverse[G0];
	Table[Tr[2K1.K2+2K1.G0.Transpose[K2].invG0],{K1,LieBasis},{K2,LieBasis}]//SparseArray];

GOMetricO[LieBasis_]:=IdentityMatrix[Length[LieBasis]]//SparseArray;

(* Natural geometry for ONB bases *)
(* This computes the natural Fubini-Study metric associated to Lie alebra elements LieBasis around the state J0. *)
(* The code works for bosons and fermions at the same time by choosing pm=+1 for bosons and pm=-1 for fermions! *)
(* ATTENTION: This code only works for J0 in standard form!!! Ortherwise, you need to use the code below GOGeometryConstOLD which is commented out*)
(* Our Lie algebra bases are constructed as ONB bases with respect to J0 in standard form. This means we can simplify the calculation temendeously in comparison to GOGeometryConstOLD.*)
(* Hoewver, we still keep the calculation. *)
GOGeometryConst[LieBasis_,J0_,pm_]:=Module[{invG0,metric,invmetric,norm2,invnorm2},
	norm2=pm Table[Tr[MatrixPower[(K.J0-J0.K),2]],{K,LieBasis}]//Chop;
	invnorm2=Table[If[x==0,0,1/x],{x,norm2}];
	metric=DiagonalMatrix[norm2]//SparseArray;
	invmetric=DiagonalMatrix[invnorm2]//SparseArray;
	Return[{metric,invmetric}]
	];

(* Natural geometry *)
(* This computes the natural Fubini-Study metric associated to Lie alebra elements LieBasis around the state J0. *)
(* The code works for bosons and fermions at the same time by choosing pm=+1 for bosons and pm=-1 for fermions! *)
(* GOGeometryConstOLD[LieBasis_,J0_,pm_]:=Module[{invG0,metric,invmetric},
	metric=pm Table[Tr[(K1.J0-J0.K1).(K2.J0-J0.K2)],{K1,LieBasis},{K2,LieBasis}]//SparseArray;
	invmetric=PseudoInverse[metric]//SparseArray;
	Return[{metric,invmetric}]
	];
*)


(* --------------------------------------------------------------------Optimization algorithm-------------------------------------------------------------------------------- *)

GOOptimize[

(* Problem-specific input arguments *)
{function_, gradfunction_}, 

(* System-specific input arguments *)
{J0_, M0_, K_,geometry_, newM_},

(* Process-specific input arguments *)
{gradtol_, Etol_, steplimit_:\[Infinity], stepcorrection_, trackall_:False}]:=

	Module[{G0, invG0, (* Initial covariance matrix *)
			Mold, Mnew, Eold, Enew,(* Updating function values *)
			grad, Normgrad, X, \[Epsilon], GenerateM, invmetric, (* Movement *)
			M0list, Elist, Normlist, diffE, diffNorm, (* Tracking values *)
			stepcount, dimM0, CorrList, order1, order2, order, keepnumber, loosenumber, donelist, stopreason,(* Tracking trajectories *) 
			resultIndex, doneE, doneElist, doneM, doneNorm, FinalE, FinalM, FinalElist, FinalNormlist}, (* Results *)
		
		dimM0=Length[M0]; CorrList=List[]; M0list=M0; donelist=List[]; Elist=List[];
		
		(* Sub-rountine for step size - this definition is always the same but depends on the stepcorrection function *)
		GenerateM[\[Epsilon]_,Mold_,Mnew_,Eold_,Enew_,X_]:=Module[{s=\[Epsilon], enew=Enew, corr=0, mnew=Mnew},
			While[enew>Eold, s=stepcorrection[s]; corr++; mnew=newM[Mold,s,X]; enew=function[mnew,J0];]; AppendTo[CorrList,corr]; Return[{mnew//SparseArray,enew}]];
		
		(* --------Iteration-------- *)
		
		(* Define function values and gradient for initial values *) 
		Mold=M0; Eold=function[#,J0]&/@Mold; {grad,Normgrad,X}=(gradfunction[#,J0,K,geometry]&/@Mold)//Transpose;
		
		(* Print[Eold]; Lucas: Why print it here?*)
		
		Elist={Eold}//Transpose; Normlist={Normgrad}//Transpose; diffE=Eold; diffNorm=Normgrad; 
		
		(* Initialise step count *)
		stepcount=0;	
		
		(* -----Main routine----- *)
		
		(* Stopping condition *)
		While[AllTrue[Normgrad,#>gradtol&] && stepcount < steplimit && AllTrue[diffE,#>Etol&],
		
			(* Choose initial step size *)
			\[Epsilon]=Normgrad/2;
			
			(* Calculate new transformation / function value *)
			Mnew=MapThread[newM,{Mold,\[Epsilon],X}]; Enew=function[#,J0]&/@Mnew;
			
			(* Sub-rountine to ensure favourable step *)
			{Mnew,Enew}=MapThread[GenerateM,{\[Epsilon],Mold,Mnew,Eold,Enew,X}]//Transpose;
		
		(* Define function values and gradient for new values *) 
		Mold=Mnew; Eold=Enew; {grad,Normgrad,X}=(gradfunction[#,J0,K,geometry]&/@Mold)//Transpose; 
		Elist=Elist//Transpose; Elist=AppendTo[Elist,Eold]//Transpose; Normlist=Normlist//Transpose; Normlist=AppendTo[Normlist,Normgrad]//Transpose;
		diffE=Abs[#[[-1]]-#[[-2]]]&/@Elist; diffNorm=Abs[#[[-1]]-#[[-2]]]&/@Normlist;

		(* Update step count *)
		stepcount++;
											
		(* Checking trajectories *)
		If[stepcount/5//IntegerQ && Length[M0list]>1 && !trackall,
			
			(* Retain most promising 20% *)
			keepnumber=If[(.2 dimM0//Round)==0,1,.2 dimM0//Round]; 
			
			order1=Ordering[Eold,All,Less]; order2=Ordering[Normgrad,All,Greater]; 
			If[keepnumber>1,order=Join[Take[order1,keepnumber],Take[order2,keepnumber]],order=Take[order1,keepnumber]];
			
			{Mold, Eold, grad, Normgrad, X, M0list, Elist, Normlist, diffE, diffNorm}=(#[[order]])&/@{Mold, Eold, grad, Normgrad, X, M0list, Elist, Normlist, diffE, diffNorm}; dimM0=Length[order];	
			];
		];
	
	(* Check overall stopping criteria *)
	If[#<gradtol,AppendTo[donelist,Position[Normgrad,#]];stopreason="Gradient norm tolerance reached"]&/@Normgrad; 
	If[stepcount >= steplimit,AppendTo[donelist,Position[Eold,Min[Eold]]];stopreason="Iteration limit reached"];
	If[#<Etol,AppendTo[donelist,Position[diffE,#]];stopreason="Function value tolerance reached"]&/@diffE;
	
	doneE=Eold[[donelist//Flatten]]; doneElist=Elist[[donelist//Flatten]]; doneM=Mold[[donelist//Flatten]]; doneNorm=Normlist[[donelist//Flatten]];
	resultIndex=Position[doneE,Min[doneE]]//Flatten;
	
	FinalE=DeleteDuplicatesBy[doneE[[resultIndex]],10^-10]; FinalM=doneM[[resultIndex]]; FinalElist=doneElist[[resultIndex]]; FinalNormlist=doneNorm[[resultIndex]];
	(* --- Output --- *)

	Return[{FinalE,FinalM,CorrList,FinalElist,Elist,FinalNormlist,stopreason}]
];


(* --------------------------------------------------------------------Purifications and the standard form-------------------------------------------------------------------- *)

(* Extracting the parameters for the standard form *)
(*GOExtractStdFormG[G_]:=Module[{NT=Dimensions[G][[1]]/2,\[CapitalOmega]T,J,RI,SWT,tra,resc,Tra,GAB,ON,TRA,check,checklist,DiagElems,Diag,Mtra,rlist}, 
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
	Return[{rlist,Mtra}]];*)
GOExtractStdFormG[G_]:=Module[{NT=Dimensions[G][[1]]/2,sqG,\[CapitalOmega]T,mat,q,t,JJ,Mtra,sign,MTra,rlist}, 
	sqG=MatrixPower[G,-1/2];
	\[CapitalOmega]T=GO\[CapitalOmega]qpqp[NT];
	mat=sqG.\[CapitalOmega]T.sqG;
	JJ=\[CapitalOmega]T.G;
	{q,t}=SchurDecomposition[mat]//Chop;
	Mtra=Transpose[sqG.q.MatrixPower[t,-1/2]];
	rlist=Table[1/2 ArcCosh[(Mtra.G.Transpose[Mtra])[[2i,2i]]],{i,NT}];
	sign=Mtra.GOinvMSp[Mtra]//Chop;
	MTra=DiagonalMatrix[Table[If[sign[[i,i]]<0&&OddQ[i],-1,1],{i,2NT}]].Mtra;
	Return[{rlist,MTra}]];
	
(*GOExtractStdForm\[CapitalOmega][\[CapitalOmega]_]:=Module[{tuples,selecttuples,values,vectors,g,invg,Mtra,rlist},
	tuples=Eigensystem[\[CapitalOmega]]//Chop//Transpose;
	selecttuples=Select[tuples,Im[#][[1]]>=0&]//Chop;
	{values,vectors}=selecttuples//Transpose;
	g=vectors.ConjugateTranspose[vectors]//Chop;
	invg=MatrixPower[Inverse[g],1/2];
	Mtra=Table[{Sqrt[2]*Re[x],Sqrt[2]*Im[x]},{x,invg.vectors}]//Flatten[#,1]&;
	rlist=1/2ArcCos[Im[#]]&/@values//N;
	Return[{rlist,Mtra}]];*)
GOExtractStdForm\[CapitalOmega][\[CapitalOmega]_]:=Module[{evalues,evectors,RealImaginaryVectors,SelectedVectors,rList,Mtra},
	{evalues,evectors}=Eigensystem[\[CapitalOmega]]//Chop;
	RealImaginaryVectors=Flatten[Table[{Re[vec],Im[vec]},{vec,evectors}],1];
	SelectedVectors=DeleteDuplicates[RealImaginaryVectors,(#1==-#2)||(#1==#2)||Norm[#2]==0&];
	Mtra=Orthogonalize[Table[vec/Norm[vec],{vec,SelectedVectors}]];
	rList=Table[1/2ArcCos[Abs[evalues[[i]]]],{i,1,Length[evalues],2}];
	Return[{rList,Mtra}]]

GOPurifyStandardGBoson[rlist_,dimp_,Gform_]:=Module[{cosh,sinh,diag,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0,n,m,Mtra},
	cosh=Flatten[Table[{Cosh[2 rr],Cosh[2 rr]},{rr,rlist}]];
	sinh=Flatten[Table[{Sinh[2 rr],-Sinh[2 rr]},{rr,rlist}]];
	
	n=Length[rlist]; m=dimp;
	
	Q14=DiagonalMatrix[Hold/@cosh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sinh]//ReleaseHold//ArrayFlatten//SparseArray;
	Q5=If[m==n,Null,IdentityMatrix[2(m-n)]//SparseArray];

	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{Q23,Q14}}]//SparseArray, G0=ArrayFlatten[{{Q14,Q23,0},{Q23,Q14,0},{0,0,Q5}}]//SparseArray];
	
	Mtra=If[Gform=="qpqp",IdentityMatrix[2(n+m)],ToExpression["GO"<>Gform<>"FROMqpqp"][n+m]];
	SparseArray[Mtra.G0.Transpose[Mtra]]];
	
GOPurifyStandardJBoson[rlist_,dimp_,Gform_]:=Module[{n,G,\[CapitalOmega]},
	G=GOPurifyStandardGBoson[rlist,dimp,Gform];
	\[CapitalOmega]=ToExpression["GO\[CapitalOmega]"<>Gform][Length[G]/2];
	SparseArray[-G.\[CapitalOmega]]];
	
GOPurifyStandard\[CapitalOmega]Fermion[rlist_,dimp_,\[CapitalOmega]form_]:=Module[{n,m,cos,sin,Q14,Q23,Q5,G0,\[CapitalOmega]0,J0,Mtra},
	cos=Table[ArrayFlatten[{{0,Cos[2 rr]},{-Cos[2 rr],0}}],{rr,rlist}];
	sin=Table[ArrayFlatten[{{0,Sin[2 rr]},{Sin[2 rr],0}}],{rr,rlist}];
	
	n=Length[rlist]; m=dimp;

	Q14=DiagonalMatrix[Hold/@cos]//ReleaseHold//ArrayFlatten//SparseArray;
	Q23=DiagonalMatrix[Hold/@sin]//ReleaseHold//ArrayFlatten//SparseArray;
	
	Q5=If[m==n,Null,ArrayFlatten[DiagonalMatrix[ConstantArray[a,m-n]]/.a -> {{0, 1},{-1, 0}}]];
	If[m==n, G0=ArrayFlatten[{{Q14,Q23},{-Q23,Q14}}]//SparseArray, 
	G0=ArrayFlatten[{{Q14,Q23,0},{-Q23,Q14,0},{0,0,Q5}}]//SparseArray];
	
	Mtra=If[\[CapitalOmega]form=="qpqp",IdentityMatrix[2(n+m)],ToExpression["GO"<>\[CapitalOmega]form<>"FROMqpqp"][(n+m)]];
	SparseArray[Mtra.G0.Transpose[Mtra]]];
	
GOPurifyStandardJFermion[rlist_,dimp_,\[CapitalOmega]form_]:=Module[{n,G,\[CapitalOmega]},
	\[CapitalOmega]=GOPurifyStandard\[CapitalOmega]Fermion[rlist,dimp,\[CapitalOmega]form];
	G=ToExpression["GOG"<>\[CapitalOmega]form][Length[\[CapitalOmega]]/2];
	SparseArray[-G.\[CapitalOmega]]];


(* --------------------------------------------------------------------------Applications------------------------------------------------------------------------------------ *)

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
	D=(M.(IdentityMatrix[n]+I J0).Transpose[M]/2)[[ResAA,ResAA]];
	logD=GOConditionalLog[D];
	Re[-Tr[D.logD]]//Chop]];
GOEoPgradBos[ResAA_]:=Function[{M,J0,K,geometry},Module[{n,nn,invM,dJ,D,dD,logDD,grad,Normgrad,X,invmetric},
	invmetric=geometry[[2]];
	n=Length[J0];
	invM=GOinvMSp[M];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,K}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logDD=GOConditionalLog[D.D];
	grad=invmetric.Table[Re[1/2 Tr[dDi.logDD]]//Chop,{dDi,dD}];
	Normgrad=Norm[grad]; X=-grad.K/Normgrad;
	{grad,Normgrad,X//SparseArray}]];
(*GOEoPgradBos[ResAA_]:=Function[{M,J0,K,geometry},Module[{n,nn,invM,dJ,D,dD,logDD,grad,Normgrad,X,metric,metricinv},
	n=Length[J0];
	invM=GOinvMSp[M];
	dJ=Table[KIi.M.J0.invM-M.J0.invM.KIi,{KIi,K}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logDD=GOConditionalLog[D.D];
	metric=GOMetricSp[K,J0];
	grad=metricinv.Table[Re[1/2 Tr[dDi.logDD]]//Chop,{dDi,dD}];
	Normgrad=Norm[grad]; X=-grad.K/Normgrad;
	{grad,Normgrad,X//SparseArray}]];*)
GOEoPgradFerm[ResAA_]:=Function[{M,J0,K,geometry},Module[{n,nn,invM,KI,dJ,D,dD,logD,grad,Normgrad,X,invmetric},
	invmetric=geometry[[2]];
	n=Length[J0]; nn=Length[J0];
	invM=Transpose[M];
	KI=Table[PadLeft[Ki,{n,n}],{Ki,K}];
	dJ=Table[M.(KIi.J0-J0.KIi).invM,{KIi,KI}];
	dD=Table[I/2 dJAAi[[ResAA,ResAA]],{dJAAi,dJ}];
	D=(M.(IdentityMatrix[n]+I J0).invM/2)[[ResAA,ResAA]];
	logD=GOConditionalLog[D];
	grad=invmetric.Table[Re[-Tr[dDi.logD]]//Chop,{dDi,dD}];
	Normgrad=Norm[grad]; X=-grad.K/Normgrad;
	{grad,Normgrad,X//SparseArray}]];
	
(* --------------- Complexity of Purification ---------------- *)
(* Functionals and gradients for bosonic and fermionic CoP *)
GOCoPBos[JT_]:=Function[{M,Jref0},Module[{invM,D},
	D=M.Jref0.GOinvMSp[M].Inverse[JT]; 
	Re[Sqrt[Total[Log[#1]^2&/@Eigenvalues[D]]/8]]]];
GOCoPFerm[JT_]:=Function[{M,Jref0},Module[{invM,D},
	D=M.Jref0.Transpose[M].Inverse[JT]; 
	Re[Sqrt[Total[I Log[#1]^2&/@Eigenvalues[D]]/8]]]];	
	
GOCoPgradBos[JT_]:=Function[{M,Jref0,K,geometry},Module[{dimA,dimB,invM,invJT,D,invD,dD,grad,Normgrad,X,invmetric},
	invmetric=geometry[[2]];
	dimB=Length[K[[1]]]; dimA=Length[M]-dimB;
	invM=GOinvMSp[M]; invJT=Inverse[JT];
	D=M.Jref0.invM.invJT; invD=Inverse[D];
	dD=Table[M.(KIi.Jref0-Jref0.KIi).invM.invJT,{KIi,K}];
	grad=invmetric.Table[Re[2Tr[GOConditionalLog[D].invD.dDi]],{dDi,dD}];
	Normgrad=Norm[grad]; X=-grad.K/Normgrad;
	{grad,Normgrad,X//SparseArray}]];
GOCoPgradFerm[JT_]:=Function[{M,Jref0,K,geometry},Module[{dimA,dimB,invM,invJT,D,invD,dD,grad,Normgrad,X,invmetric},
	invmetric=geometry[[2]];
	dimB=Length[K[[1]]]; dimA=Length[M]-dimB;
	invM=Transpose[M]; invJT=Inverse[JT];
	D=M.Jref0.invM.invJT; invD=Inverse[D];
	dD=Table[M.(KIi.Jref0-Jref0.KIi).invM.invJT,{KIi,K}];
	grad=invmetric.Table[Re[-2Tr[GOConditionalLog[D].invD.dDi]],{dDi,dD}];
	Normgrad=Norm[grad]; X=-grad.K/Normgrad;
	{grad,Normgrad,X//SparseArray}]];


End[]

EndPackage[]
