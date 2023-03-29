(* ::Package:: *)

(* ::Section:: *)
(*Info*)


(* :Author: Kirill Belov *)


(* ::Section:: *)
(*Begin package*)


BeginPackage["WLGPNTeam`GeologyTools`DatasetExtensions`"]


(* ::Section:: *)
(*Names*)


AddColumn::usage = 
"AddColumn[dataset, dataset2, keys]"


(* ::Section:: *)
(*Begin private*)


Begin["`Private`"]


(* ::Section:: *)
(*AddColumn*)


AddColumn /: Dot[dataset: {__Association}, AddColumn[key_, column_, dataset2_]] := 
	Map[SelectFirst[dataset2, #]&] @ dataset[[All, key]]


(* ::Section:: *)
(*End private*)


End[] (*`Private`*)


(* ::Section:: *)
(*End package*)


EndPackage[] (*GPNTools`DatasetExtensions`*)
