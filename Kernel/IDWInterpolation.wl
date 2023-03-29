(* ::Package:: *)

(* ::Section:: *)
(*Info*)


(* :Package: GPNTools *)
(* :Context: GPNTools`IDWInterpolation` *)
(* Version: 0.1.1 *)
(* Developer: Kirill Belov *)


(* ::Section:: *)
(*Begin package*)


BeginPackage["WLGPNTeam`GeologyTools`IDWInterpolation`"]


(* ::Section:: *)
(*Names declaration*)


ClearAll[IDWInterpolation]


IDWInterpolation::usage = 
"IDWInterpolation[array3D, opts]
returns the special function IDWInterpolatedFunction[..]
Options: {Neighbors -> Integer, Radius -> ?NumericQ, Delta -> ?NumericQ, Beta -> ?NumericQ}"


ClearAll[IDWInterpolatedFunction]


IDWInterpolatedFunction::usage = 
"IDWInterpolatedFunction[{data}][x, y]"


(* ::Section:: *)
(*Begin private*)


Begin["`Private`"]


(* ::Section:: *)
(*Internal functions*)


(* ::Text:: *)
(*Compiled function - fast calc of the IDW*)


$IDWFunction = Compile[{
	{points, _Real, 2}, 
	{radius}, {delta}, {beta}, 
	{x}, {y}, {n, _Integer}
}, 
	Block[{
		values, 
		weights, 
		dists, 
		mdist
	}, 
		values = points[[All, 3]]; 
		dists = Table[
			(x - p[[1]])^2 + (y - p[[2]])^2, 
			{p, points}
		]; 
		mdist = Sort[dists][[If[n <= Length[points], n, Length[points]]]];
		weights = Table[
			If[dists[[i]] == 0, Return[values[[i]]]]; 
			If[dists[[i]] <= mdist, 1 / (dists[[i]] + delta^2)^beta, 0.0], 
			{i, 1, Length[dists]}
		]; 
		Return[Total[values * weights] / Total[weights]]
	], 
	CompilationOptions -> {"ExpressionOptimization" -> True}, 
	Parallelization -> True, 
	RuntimeAttributes -> {Listable}
];


(* ::Text:: *)
(*Finding the block for the point (x, y)*)


$IDWGetBlock[{{xmin_?NumericQ, xmax_?NumericQ}, 
	{ymin_?NumericQ, ymax_?NumericQ}}, radius_?NumericQ][x_?NumericQ, y_?NumericQ] := 
With[{xmean = Mean[{xmin, xmax}], ymean = Mean[{ymin, ymax}]}, 
	{Floor[(x - xmean)/radius], Floor[(y - ymean)/radius]}
]


(* ::Text:: *)
(*Finding the main block and nears*)


$IDWGetBlocks[{{xmin_?NumericQ, xmax_?NumericQ}, {ymin_?NumericQ, ymax_?NumericQ}}, 
	radius_?NumericQ][x_?NumericQ, y_?NumericQ] := 
With[{xmean = Mean[{xmin, xmax}], ymean = Mean[{ymin, ymax}]}, 
	{
		{Floor[(x - xmean - radius/2)/radius], Floor[(y - ymean - radius/2)/radius]}, 
		{Floor[(x - xmean - radius/2)/radius], Floor[(y - ymean + radius/2)/radius]}, 
		{Floor[(x - xmean + radius/2)/radius], Floor[(y - ymean - radius/2)/radius]}, 
		{Floor[(x - xmean + radius/2)/radius], Floor[(y - ymean + radius/2)/radius]}
	}
]


(* ::Text:: *)
(**)


(* ::Section:: *)
(*Main functions*)


(* ::Text:: *)
(*IDWInterpolation*)


(* ::Text:: *)
(*Options*)


Options[IDWInterpolation] = {
	"Neighbors" -> 6, 
	"Radius" -> Automatic, 
	"Delta" -> Automatic, 
	"Beta" -> 2
};


(* ::Text:: *)
(*Main defenition*)


IDWInterpolation[mesh: {{_?NumericQ, _?NumericQ, _?NumericQ}..}, 
	radius_?NumericQ, delta_?NumericQ, beta_?NumericQ, n_Integer] := 
Block[{
	meshes = <||>, 
	range = {MinMax[mesh[[All, 1]]], MinMax[mesh[[All, 2]]]}, 
	blocks
}, 
	Table[
		blocks = $IDWGetBlocks[range, radius][p[[1]], p[[2]]]; 
		Table[
			If[Not[KeyExistsQ[meshes, block]], meshes[block] = {}];
			meshes[block] = Append[meshes[block], p], 
			{block, blocks}
		]; , 
		{p, mesh}
	];
	IDWInterpolatedFunction[{meshes, range, radius, delta, beta, n}]
]


(* ::Text:: *)
(*Advanced defenition with opts*)


IDWInterpolation[mesh: {{_?NumericQ, _?NumericQ, _?NumericQ}..}, OptionsPattern[]] := 
Block[{
	nneighbors = OptionValue["Neighbors"], 
	radius = OptionValue["Radius"], 
	delta = OptionValue["Delta"], 
	beta = OptionValue["Beta"]
}, 
	Which[
		NumericQ[radius], Null, 
		radius === Automatic && IntegerQ[nneighbors], 
			radius = 2 * Sqrt[nneighbors / 
				(Length[mesh] / ((Max[mesh[[All, 1]]] - Min[mesh[[All, 1]]]) * 
				(Max[mesh[[All, 2]]] - Min[mesh[[All, 2]]])))]
	]; 
	Which[
		NumericQ[delta], Null, 
		delta === Automatic, 
			delta = (Max[mesh[[All, 3]]] - Min[mesh[[All, 3]]]) / 1000.0
	]; 
	IDWInterpolation[mesh, radius, delta, beta, nneighbors]
]


(* ::Text:: *)
(*IDWInterpolatedFunction*)


(* ::Text:: *)
(*As function of two arguments*)


IDWInterpolatedFunction[{
	meshes_Association?AssociationQ, 
	range: {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}, 
	radius_?NumericQ, 
	delta_?NumericQ, 
	beta_?NumericQ, 
	n_Integer
}][x_?NumericQ, y_?NumericQ] := 
With[{
	
}, 
	$IDWFunction[meshes[$IDWGetBlock[range, radius][x, y]], 
		radius, delta, beta, x, y, n
	]
]


(* ::Text:: *)
(*As function of the one argument - list of two coordinates*)


IDWInterpolatedFunction[{
	meshes_Association?AssociationQ, 
	range: {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}, 
	radius_?NumericQ, 
	delta_?NumericQ, 
	beta_?NumericQ, 
	n_Integer
}][{x_?NumericQ, y_?NumericQ}] := 
With[{
	
}, 
	$IDWFunction[meshes[$IDWGetBlock[range, radius][x, y]], 
		radius, delta, beta, x, y, n
	]
]


(* ::Text:: *)
(*As listable function*)


IDWInterpolatedFunction[{
	meshes_Association?AssociationQ, 
	range: {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}, 
	radius_?NumericQ, 
	delta_?NumericQ, 
	beta_?NumericQ, 
	n_Integer
}][points: {{_?NumericQ, _?NumericQ}..}] := 
With[{
	
}, 
	$IDWFunction[
		Table[meshes[$IDWGetBlock[range, radius][p[[1]], p[[2]]]], {p, points}], 
		radius, delta, beta, points[[All, 1]], points[[All, 2]], n
	]
]


(* ::Text:: *)
(*Params*)


(function: IDWInterpolatedFunction[{__}])[param: 
	"Range" | "Radius" | "Delta" | "Beta" | "Neighbors"
] := Switch[param, 
	"Range", function[[1, 2]], 
	"Radius", function[[1, 3]], 
	"Delta", function[[1, 4]], 
	"Beta", function[[1, 5]], 
	"Neighbors", function[[1, 6]]
] 


(* ::Text:: *)
(*Function representation*)


IDWInterpolatedFunction /: 
MakeBoxes[function: IDWInterpolatedFunction[{
	meshes_Association?AssociationQ, 
	range: {{_?NumericQ, _?NumericQ}, {_?NumericQ, _?NumericQ}}, 
	radius_?NumericQ, 
	delta_?NumericQ, 
	beta_?NumericQ, 
	n_Integer
}], form: (StandardForm | OutputForm)] := 
Module[{above, below, icon = Null}, 
	icon = Magnify[ContourPlot[
		function[x, y], 
		{x, range[[1, 1]], range[[1, 2]]}, 
		{y, range[[2, 1]], range[[2, 2]]}, 
		PlotPoints -> 3, 
		Frame -> False, 
		Axes -> False, 
		ImageSize -> Tiny, 
		PlotRange -> All, 
		AspectRatio -> Full, 
		PlotTheme -> "Monochrome"
	], 1/2]; 
	
	above = {
		{
			BoxForm`SummaryItem[{
				"Range: ", 
				Row[{
					"x: ", NumberForm[range[[1]]], 
					" y: ", NumberForm[range[[2]]]
				}]
			}], 
			SpanFromLeft
		}, 
		{
			BoxForm`SummaryItem[{"Radius: ", radius}], 
			BoxForm`SummaryItem[{"Delta: ", delta}], 
			BoxForm`SummaryItem[{"Beta: ", beta}], 
			BoxForm`SummaryItem[{"Neighbors: ", n}], 
			SpanFromLeft
		}
	};
	
	below = {
		
	}; 
	
    BoxForm`ArrangeSummaryBox[
		IDWInterpolatedFunction, 
		function, 
		icon, 
		above, 
		below, 
		form, 
		"Interpretable" -> Automatic
	]
]


(* ::Section:: *)
(*End private*)


End[] (*`Private`*)


(* ::Section:: *)
(*End package*)


EndPackage[] (*GPNTools`IDWInterpolation`*)
