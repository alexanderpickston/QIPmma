(* ::Package:: *)

(* Wolfram Language Package *)

(* Created by Alexander Pickston, adapted from code originally written by Massimiliano Proietti *)

(* alexpickston@gmail.com *)

BeginPackage["gsl`"]

CustomGraph::usage="CustomGraph[edges_]
From a list of edges, create a graph using. To make a graph equivalent to a GHZ state, then edges={{1,2},{2,3}}"

(* graph styling *)
graphstyle = {
   	VertexSize -> {0.4},
   	VertexLabels -> {Placed["Name", {1/2*1.05, 1/2*1.05}]},
   	ImageSize -> {146.99479166666538, Automatic},
   	VertexStyle -> {Directive[EdgeForm[{Thick, Opacity[1],Blue}], Blue]},
   	VertexLabelStyle -> Directive[White, FontFamily -> "IBM Plex Mono", 20],
   	EdgeStyle -> Directive[Black, Thick, Opacity[1]],
   	GraphLayout -> "CircularEmbedding"
    };

(* folded Kronecker product *)
Kron := KroneckerProduct[##]&

(* basis states *)
h = {{1}, {0}};
v = {{0}, {1}};
d = 1/Sqrt[2] (h + v);
a = 1/Sqrt[2] (h - v);
r = 1/Sqrt[2] (h + I*v);
l = 1/Sqrt[2] (h - I*v);

(* Partial Trace - adapted from Mark S. Tame - http://library.wolfram.com/infocenter/MathSource/5571/ *)
SwapParts[expr_, pos1_, pos2_] := 
    
    ReplacePart[#, #, {pos1, pos2}, {pos2, pos1}] &[expr];
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
            Take[M[[p + 1, All]], {2, 2^n, 2}]];
            (* Pick row p/p+1 and take elements 1/2 through 2^n in steps of 2 - sum thoese and append to zero matrix - do for all rows *)
            ],

        For[j = 0, j < (n - k), j++,(* if not - permute accordingly *) b = {0};
        For[i = 1, i < 2^n + 1, i++,
        If[(Mod[(IntegerDigits[i - 1, 2, n][[n]] + 
            IntegerDigits[i - 1, 2, n][[n - j - 1]]), 2]) == 1 && 
            Count[b, i]  == 0, 
            Permut = {i, (FromDigits[
                SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 1}], 2] + 1)};
                b = Append[b, (FromDigits[
                SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 1}], 2] + 1)];
                c = Range[2^n];
                perm = 
                SwapParts[c, {i}, {(FromDigits[SwapParts[(IntegerDigits[i - 1, 2, n]), {n}, {n - j - 1}], 2] + 1)}];
                M = M[[perm, perm]];]];
        
    (* and now trace out last system *)
    TrkM = {};
    For[p = 1, p < 2^n + 1, p = p + 2,
    TrkM = Append[TrkM, 
        Take[M[[p, All]], {1, 2^n, 2}] + 
        Take[M[[p + 1, All]], {2, 2^n, 2}]];
        ]
        ]
        ];
    ];
    
    Clear[Qubits, n, M, k, Permut, perm, b, c]; 
    TrkM];

Module[{input, out},
    CustomGraph[edges_] := (
    
        (* edges should be a string*)
        (* should be of the form: 
        edges={{1,2},{2,3}} *)

        input = Table[edges[[i, 1]] \[UndirectedEdge] edges[[i, 2]]
          , {i, Length@edges}];

        out = Graph[input, graphstyle];
    Return[out])
];

Module[{subGraph,diffGraph,complementGraph,out,g,vertexCoordinates},
	LCQubit[graph_,vertex_]:=(

		(* get the vertex coordinates of the original graph *)
		vertexCoordinates=GraphEmbedding[graph];

		(* apply the vertex coordinates to the original graph *)
		g=Graph[VertexList[graph],EdgeList[graph],VertexCoordinates->vertexCoordinates];

		(* select the subgraph generated by the vertex and its neighbours *)
		subGraph=Subgraph[graph,AdjacencyList[graph,vertex]];

		(* complement the sub graph *)
		complementGraph=GraphComplement[subGraph];

		(* remove the starting subgraph *)
		diffGraph=GraphDifference[graph,subGraph];

		(* union the new subGraph with the remaining main graph *)
		out=GraphUnion[diffGraph,complementGraph];
		out=Graph[out,graphstyle,VertexCoordinates->vertexCoordinates];

	(* return the LC-equivalent graph *)
	Return[out];);
];

Module[{perm,prmList,noDuplicates,g,vertexCoordinates,out},
	LCOrbit[graph_]:=(

		(* get the vertex coordinates of the original graph *)
		vertexCoordinates=GraphEmbedding[graph];

		(* apply the vertex coordinates to the original graph *)
		g=Graph[VertexList[graph],EdgeList[graph],VertexCoordinates->vertexCoordinates];

		(* module operations *)
		perm=Permutations[Range[VertexCount[graph]],VertexCount[graph]];
		prmList=FoldList[LCQubit,graph,#]&/@perm//Flatten;
		noDuplicates=DeleteDuplicates[prmList,IsomorphicGraphQ];

		(* ensure styling is correct *)
		out=Graph[noDuplicates,graphstyle,VertexCoordinates->vertexCoordinates];

	(* return the LC-equivalent graph orbit *)
	Return[noDuplicates];)
];

Module[{vertexCoordinates,g,edgeDelete,edgeDeleteGraph,vertexList,vertexDeleted,ordering,out},
	Zmeasurement[graph_,vertex_]:=(

		(* get the vertex coordinates of the original graph *)
		vertexCoordinates=GraphEmbedding[graph];

		(* apply the vertex coordinates to the original graph *)
		g=Graph[VertexList[graph],EdgeList[graph],graphstyle,VertexCoordinates->vertexCoordinates];

		(* module operation *)
		(* finding the complement between all edges and those edges which join to the vertex specified in the function *)
		edgeDelete=Complement[EdgeList[g],EdgeList[g,vertex\[UndirectedEdge]_]];
		edgeDeleteGraph=Graph[edgeDelete];
		vertexList=VertexList[edgeDeleteGraph];

		(* need to work out what edges have been deleted as some vertices will also be removed *)
		(* only showing vertices that still posses an edge *)
		vertexDeleted=Complement[Range@Length@vertexCoordinates,vertexList];
		(* correct re-formatting of variable to parse into the Delete[] function *)
		vertexDeleted={#}&/@vertexDeleted;
		(* new vertex co-ordinates with deleted vertices removed *) 
		vertexCoordinates=Delete[vertexCoordinates,vertexDeleted];

		(* need to get the order correct *)
		(* rearrange the remaining vertexCoordinates with the VertexList[] for the current graph *)
		ordering=Ordering@VertexList[edgeDelete];
		vertexCoordinates=vertexCoordinates[[#]]&/@ordering;

		(* ensure styling is correct *)
		out=Graph[edgeDelete,graphstyle,VertexCoordinates->vertexCoordinates];

	Return[out];);
];

Module[{vertexCoordinates,g,edgeDelete,edgeDeleteGraph,vertexList,vertexDeleted,ordering,out},
	ZmeasurementOLD[graph_,vertex_]:=(

		(* get the vertex coordinates of the original graph *)
		vertexCoordinates=GraphEmbedding[graph];

		(* apply the vertex coordinates to the original graph *)
		g=Graph[VertexList[graph],EdgeList[graph],graphstyle,VertexCoordinates->vertexCoordinates];

		(* module operation *)
		(* finding the complement between all edges and those edges which join to the vertex specified in the function *)
		edgeDelete=Complement[EdgeList[g],EdgeList[g,vertex\[UndirectedEdge]_]];
		edgeDeleteGraph=Graph[edgeDelete];
		vertexList=VertexList[edgeDeleteGraph];

		(* need to work out what edges have been deleted as some vertices will also be removed *)
		(* only showing vertices that still posses an edge *)
		vertexDeleted=Complement[Range@Length@vertexCoordinates,vertexList];
		(* correct re-formatting of variable to parse into the Delete[] function *)
		vertexDeleted={#}&/@vertexDeleted;
		(* new vertex co-ordinates with deleted vertices removed *) 
		vertexCoordinates=Delete[vertexCoordinates,vertexDeleted];

		(* need to get the order correct *)
		(* rearrange the remaining vertexCoordinates with the VertexList[] for the current graph *)
		ordering=Ordering@VertexList[edgeDelete];
		vertexCoordinates=vertexCoordinates[[#]]&/@ordering;

		(* ensure styling is correct *)
		out=Graph[edgeDelete,graphstyle,VertexCoordinates->vertexCoordinates];

	Return[out];);
];