(* ::Package:: *)

(* ::Title:: *)
(*QiPmma*)


(* Wolfram Language Package *)

(* Created by Alexander Pickston *)

(* alexpickston@gmail.com *)


BeginPackage["QIPmma`"]


(* docs *)
Fidelity::usage="Fidelity[\!\(\*
StyleBox[\"rho\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"target\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"rho\",\nFontColor->RGBColor[0, 0, 1]]\)]
Computes the fidelity between a target density operator and a density operator defined by the user";


Fidelity::usage="FidelityPure[\!\(\*
StyleBox[\"state\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"target\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"state\",\nFontColor->RGBColor[0, 0, 1]]\)]
Computes the fidelity between a target pure quantum state and a pure quantum state defined by the user";


Purity::usage="Purity[\!\(\*
StyleBox[\"rho\",\nFontColor->RGBColor[0, 0, 1]]\)]
Calculate the purity of a quantum state rho"


DensityMatrix::usage="DensityMatrix[\!\(\*
StyleBox[\"state\",\nFontColor->RGBColor[0, 0, 1]]\)]
Construct the density matrix of a given state"


Kron::usage="Kron[\!\(\*
StyleBox[\"#1\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"#2\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"#3\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"etc\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"...\",\nFontColor->RGBColor[0, 0, 1]]\)]
Calculate the Kronecker product of more than two vectors/matrices"


StateMeasurementPure::usage="StateMeasurementPure[\!\(\*
StyleBox[\"state\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"operator\",\nFontColor->RGBColor[0, 0, 1]]\)]
Calculates the expectation value of an operator for a given state"


StateMeasurement::usage="StateMeasurement[\!\(\*
StyleBox[\"rho\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\",\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
StyleBox[\"operator\",\nFontColor->RGBColor[0, 0, 1]]\)]
Calculates measurement outcome using the Born rule"


EigenSolve::usage="EigenSolve[\!\(\*
StyleBox[\"matrix\",\nFontColor->RGBColor[0, 0, 1]]\)]
Returns the eigenvalues and eigenvectors of the given matrix"


GetProjectionAngles::usage="GetProjectionAngles[\!\(\*
StyleBox[\"projector\",\nFontColor->RGBColor[0, 0, 1]]\)]
Given a vector (or projector) which is a projection measurement in a certain basis, return the required measurement angles for a measurement stage containing a QWP \[Rule] HWP \[Rule] PBS"


Begin["`Private`"]


(* ::Subsection:: *)
(*Variables*)


(* ::Subsubsection::Closed:: *)
(*Basis states*)


h = {{1}, {0}};
v = {{0}, {1}};
d = 1/Sqrt[2] (h + v);
a = 1/Sqrt[2] (h - v);
r = 1/Sqrt[2] (h + I*v);
l = 1/Sqrt[2] (h - I*v);

hh = Kron[h, h]; hv = Kron[h, v]; hd = Kron[h, d];
ha = Kron[h, a]; hr = Kron[h, r]; hl = Kron[h, l];
vh = Kron[v, h]; vv = Kron[v, v]; vd = Kron[v, d];
va = Kron[v, a]; vr = Kron[v, r]; vl = Kron[v, l];
dh = Kron[d, h]; dv = Kron[d, v]; dd = Kron[d, d];
da = Kron[d, a]; dr = Kron[d, r]; dl = Kron[d, l];
ah = Kron[a, h]; av = Kron[a, v]; ad = Kron[a, v];
aa = Kron[a, a]; ar = Kron[a, r]; al = Kron[a, l];
rh = Kron[r, h]; rv = Kron[r, v]; rd = Kron[r, v];
ra = Kron[r, a]; rr = Kron[r, r]; rl = Kron[r, l];
lh = Kron[l, h]; lv = Kron[l, v]; ld = Kron[l, v];
la = Kron[l, a]; lr = Kron[l, r]; ll = Kron[l, l];


(* ::Subsubsection:: *)
(*Single qubit operations*)


s0 = {{1, 0}, {0, 1}};
sx = {{0, 1}, {1, 0}};
sy = {{0, -I}, {I, 0}};
sz = {{1, 0}, {0, -1}}; 
Hada = 1/Sqrt[2] {{1, 1}, {1, -1}};
LCNeighborhood = Exp[(s0*\[Pi])/4] {{1, 0}, {0, -s0}};
LCVertex = MatrixExp[(-s0*\[Pi]/4)*(sx)];


(* ::Subsubsection::Closed:: *)
(*Controlled gates*)


cNot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cH = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1/\[Sqrt]2, 1/\[Sqrt]2}, {0, 0, 1/\[Sqrt]2, -(1/\[Sqrt]2)}};
cX = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cY = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -I}, {0, 0, I, 0}};
cZ = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
cI = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
cSwap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};


(* ::Subsection:: *)
(*Functions*)


Fidelity[dm_, idealdm_] := Tr[MatrixPower[MatrixPower[idealdm, 1/2].dm.MatrixPower[idealdm, 1/2], 1/2]]^2 // Chop;


FidelityPure[x_, y_] := Sqrt[x\[ConjugateTranspose].y.x];


Purity[\[Rho]_]:=Tr[\[Rho].\[Rho]];


DensityMatrix:=#.#\[ConjugateTranspose]&;


Kron:=KroneckerProduct[##]&;


StateMeasurementPure[state_,operator_]:=ConjugateTranspose@state.operator.state;


StateMeasurement[rho_,operator_]:=Tr[operator.rho];


ProjectionMeasurement[x_, y_] := 
  FullSimplify[Abs[x\[ConjugateTranspose].y]^2][[1, 1]];


GetAngle[x_]:=N@180*FullSimplify@ArcTan[x[[2]]/x[[1]]]/\[Pi];


EigenSolve[matrix_] := 
	Block[{eigenvalues, eigenvectors},
	{eigenvalues, eigenvectors}={#[[1]],#[[2]]}&/@Eigensystem[matrix];
	{eigenvalues, Normalize /@ eigenvectors}
];


GetProjectionAngles[proj_]:=
	Module[{\[Theta],out},

		\[Theta]=Chop@NMinimize[{Abs[(HWP[x].QWP[y].proj)[[2, 1]]],-\[Pi]/4<=x<=\[Pi]/4},{x, y}][[2,{1,2},2]];
		out=Chop@FullSimplify@((\[Theta]/\[Pi])*180);
		
    Return[out]
];


GHZ[nQubit_]:=1/Sqrt[2](Kron@@ConstantArray[h,nQubit]+Kron@@ConstantArray[v,nQubit])/;nQubit>=2


BitFlip[p_, \[Rho]_] := p.cX.\[Rho].cX + (1 - p) \[Rho];


PhaseFlip[p_, \[Rho]_] := p*sz.\[Rho].sz + (1 - p) \[Rho];


BitPhaseFlip[p_, \[Rho]_] := p*sy.\[Rho].sy + (1 - p) \[Rho];


DepolControl[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(cX.\[Rho].cX + cY.\[Rho].cY + cZ.\[Rho].cZ);


Depol[\[Eta]_, \[Rho]_] := (1 - \[Eta]) \[Rho] + \[Eta] s0/2;


Depol1[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(sx.\[Rho].sx + sy.\[Rho].sy + sz.\[Rho].sz);


Depol2[p_, \[Rho]_] := p/2*cI + (1 - p) \[Rho];


DepolMemory[\[Eta]_, \[Rho]_] := \[Eta] \[Rho] + (1 - \[Eta]) s0/2;


QWP[t_]:=1/Sqrt[2] ({
 {1+I Cos[2t], I Sin[2t]},
 {I Sin[2t], 1-I Cos[2t]}
});


HWP[t_]:=Exp[(I \[Pi])/2]({
 {Cos[2t], Sin[2t]},
 {Sin[2t], -Cos[2t]}
});


End[]


EndPackage[]


chars = Characters@"QIPmma library has been loaded successfully. Have fun!";
list = Table[
   Style[chars[[i]], 
    Blend[{Green,Blue}, (i - 1)/(Length@chars - 1)], 
    Bold], {i, 1, Length@chars}];
Apply[Print, list]
