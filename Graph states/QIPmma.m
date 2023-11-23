(* __Functions__ *)

Fidelity[dm_, idealdm_] := Tr[MatrixPower[MatrixPower[idealdm, 1/2].dm.MatrixPower[idealdm, 1/2], 1/2]]^2 // Chop;

FidelityPure[x_, y_] := Sqrt[x\[ConjugateTranspose].y.x];
Purity[\[Rho]_]:=Tr[\[Rho].\[Rho]];

DensityMatrix:=#.#\[ConjugateTranspose]&
Kron:=KroneckerProduct[##]&
Kronk = Fold[KroneckerProduct];

neighborhoodLC = Exp[(I*\[Pi])/4] {{1, 0}, {0, -I}}
vertexLC = MatrixExp[(-I*\[Pi]/4)*(sx)]

StateMeasurementPure[state_,operator_]:=ConjugateTranspose@state.operator.state
StateMeasurement[rho_,operator_]:=Tr[operator.rho]

ProjectionMeasurement[x_, y_] := 
  FullSimplify[Abs[x\[ConjugateTranspose].y]^2][[1, 1]];

GetAngle[x_] := N@180*FullSimplify@ArcTan[x[[2]]/x[[1]]]/\[Pi];

EigenSolve[matrix_] := Block[{eigenvalues, eigenvectors},
  {eigenvalues, eigenvectors} = {#[[1]], #[[2]]} & /@ 
    Eigensystem[matrix];
  {eigenvalues, Normalize /@ eigenvectors}]

  GetProjectionAngles[
    proj_] :=
   {\[Theta] = 
     Chop@NMinimize[{Abs[(HWP[x].QWP[y].proj)[[2, 1]]], -\[Pi]/4 <= 
          x <= \[Pi]/4}, {x, y}][[2, {1, 2}, 2]],
    Chop@FullSimplify@((\[Theta]/\[Pi])*180)}

(*__BasisStates__*)

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

s0 = {{1, 0}, {0, 1}};
sx = {{0, 1}, {1, 0}};
sy = {{0, -I}, {I, 0}};
sz = {{1, 0}, {0, -1}}; 
Hada = 1/Sqrt[2] {{1, 1}, {1, -1}};

GHZ[nQubit_]:=1/Sqrt[2] (Kron@@ConstantArray[h,nQubit]+Kron@@ConstantArray[v,nQubit])/;nQubit>=2

(* __Gates__ *)

cNot = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cH = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1/\[Sqrt]2, 1/\[Sqrt]2}, {0, 0, 1/\[Sqrt]2, -(1/\[Sqrt]2)}};
cX = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 0}};
cY = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, -I}, {0, 0, I, 0}};
cZ = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, -1}};
cI = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
swap = {{1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}};

(* __MoreFucntions__ *)
 
BitFlip[p_, \[Rho]_] := p.cX.\[Rho].cX + (1 - p) \[Rho];
PhaseFlip[p_, \[Rho]_] := p*sz.\[Rho].sz + (1 - p) \[Rho];
BitPhaseFlip[p_, \[Rho]_] := p*sy.\[Rho].sy + (1 - p) \[Rho];
 
DepolControl[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(cX.\[Rho].cX + cY.\[Rho].cY + cZ.\[Rho].cZ);
 
Depol[\[Eta]_, \[Rho]_] := (1 - \[Eta]) \[Rho] + \[Eta] idn/2;
 
Depol1[p_, \[Rho]_] := (1 - 3/4*p) \[Rho] + p/4*(sx.\[Rho].sx + sy.\[Rho].Y + sz.\[Rho].sz);
 
Depol2[p_, \[Rho]_] := p/2*cI + (1 - p) \[Rho];
 
(* Dephase[{\[Alpha]x_, \[Alpha]y_, \[Alpha]z_}, p_, \[Rho]_] := (1 - p/2) \[Rho] + p/2 (\[Alpha]x sx.\[Rho].sx + \[Alpha]y sy.\[Rho].sy + \[Alpha]z \sz.\[Rho].sz); *)
 
DepolMemory[\[Eta]_, \[Rho]_] := \[Eta] \[Rho] + (1 - \[Eta]) idn/2;

QWP[t_]:=1/Sqrt[2] ({
 {1+I Cos[2t], I Sin[2t]},
 {I Sin[2t], 1-I Cos[2t]}
});

HWP[t_]:=Exp[(I \[Pi])/2]({
 {Cos[2t], Sin[2t]},
 {Sin[2t], -Cos[2t]}
});

TraceSystem[dm_, sys_] := 
 Block[{Qubits, TrkM, n, M, k, Permut, perm, b, c, p},
  Qubits = Reverse[Sort[sys]];(* 
  rearrange the list of qubits to be traced out *)
  TrkM = dm;
  
  (* For all qubits to be traced out... *)
  
  For[q = 1, q <= Dimensions[Qubits][[1]], q++,
   n = Log[2, (Dimensions[TrkM][[1]])]; (* 
   dimensions of original system *)
   M = TrkM;
   k = Qubits[[q]];
   
   If[k == n,(* if tracing the last system *)
    TrkM = {};
    For[p = 1, p < 2^n + 1, p = p + 2,
     TrkM = 
      Append[TrkM, 
       Take[M[[p, All]], {1, 2^n, 2}] + 
        Take[M[[p + 1, All]], {2, 2^n, 2}]];(* Pick row p/
     p+1 and take elements 1/2 through 2^n in steps of 2 - 
     sum thoese and append to zero matrix - 
     do for all rows *)
      ],
    
    For[j = 0, j < (n - k), j++,(* if not - permute accordingly *)
 
         b = {0};
     For[i = 1, i < 2^n + 1, i++,
      If[(Mod[(IntegerDigits[i - 1, 2, n][[n]] + 
             IntegerDigits[i - 1, 2, n][[n - j - 1]]), 2]) == 1 && 
        Count[b, i]  == 0, 
       Permut = {i, (FromDigits[
            SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 
               1}], 2] + 1)};
       b = 
        Append[b, (FromDigits[
            SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 
               1}], 2] + 1)];
       c = Range[2^n];
       perm = 
        SwapParts[
         c, {i}, {(FromDigits[
             SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 
                1}], 2] + 1)}];
       
       M = M[[perm, perm]];
       
        ]    
      ];
     (* and now trace out last system *)
     TrkM = {};
     For[p = 1, p < 2^n + 1, p = p + 2,
      TrkM = 
        Append[TrkM, 
         Take[M[[p, All]], {1, 2^n, 2}] + 
          Take[M[[p + 1, All]], {2, 2^n, 2}]];
      ]
        ]
    ];
   
   ];
  Clear[Qubits, n, M, k, Permut, perm, b, c];
  TrkM
  ];

Coding[n_] :=
 
 "|" <> # <> "\[RightAngleBracket]" & /@ 
  Table[IntegerString[i - 1, 2, n], {i, 1, 2^n}];