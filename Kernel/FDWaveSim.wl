(* ::Package:: *)

(* ::Title:: *)
(*FDWaveSim*)


(* ::Section:: *)
(*Info*)


(* :Title: FDWaveSim *)
(* :Context: GPNTools`FDWaveSim` *)
(* :Version: 0.1.0 *)
(* :Keywords: Waves; Finite Difference Methods; Compilation *)
(* :Authors: Anton Ekimenko; Kirill Belov *)
(* :Developer: Kirill Belov *)
(* :Email: KirillBelovTest@gmail.com *)
(* :Description: 
	
*)


(* ::Section:: *)
(*Begin package*)


BeginPackage["WLGPNTeam`GeologyTools`FDWaveSim`"]


(* ::Section:: *)
(*Unprotect and clear*)


Unprotect["`*"]
ClearAll["`*"]


(* ::Section:: *)
(*Public names declaration*)


LWaveSimulation::usage = 
"LWaveSimulation[lv, rho, signal, {xs, ys}, options]"


(* ::Section:: *)
(*Begin private*)


Begin["`Private`"]


(* ::Section:: *)
(*Compiled functions*)


LWaveSimulationCompiled[target: "MVM" | "C"] := LWaveSimulationCompiled[target] = 
	Compile[{
		{lv, _Real, 2}, 
		{rho, _Real, 2},
		{signal, _Real, 1},  
		{xs, _Integer}, 
		{ys, _Integer}
	}, 
		Module[{
			vx = 0.0 * lv, 
			vy = 0.0 * lv, 
			p = 0.0 * lv, 
			l = rho * lv * lv, 
			px = 0.0, py = 0.0, 
			dx, dy, dt, 
			nx, ny, 
			vxx = 0.0, vyy = 0.0, 
			c1, c2, c3, c4, 
			f0 = 5.0, ngpConst = 10, clfNum = 0.5
		}, 
			(*get dimension of the target environment*)
			nx = Length[p[[1]]]; 
			ny = Length[p]; 

			dx = Min[lv] / (2 * f0 * ngpConst);
			dy = dx;
			dt = dx / Max[lv] * clfNum;

			(*numeric method consts*)
			c1 = 9 / (8*dx);
			c2 = 1 / (24*dx);
			c3 = 9 / (8*dy);
			c4 = 1 / (24*dy);
	
			Table[
				Do[
					Do[
						(*Calculating spatial derivative*)
						px = c1 * (p[[ky, kx + 1]] - p[[ky, kx]]) - 
							c2 * (p[[ky, kx + 2]] - p[[ky, kx - 1]]);
					
						py = c3 * (p[[ky + 1, kx]] - p[[ky, kx]]) - 
							c4 * (p[[ky + 2, kx]] - p[[ky - 1, kx]]);
				
						(* Update velocity*)
						vx[[ky, kx]] = vx[[ky, kx]] - dt/rho[[ky, kx]] * px;
						vy[[ky, kx]] = vy[[ky, kx]] - dt/rho[[ky, kx]] * py;,
				
						{ky, 5, ny - 5, 1}
					];, 
					{kx, 5, nx - 5, 1}
				];
		
				(*Inject source wavelet*)
				p[[ys, xs]] = p[[ys, xs]] + q;
		
				(*Update pressure*)
				Do[
					Do[
						(*Calculating spatial derivative*)
						vxx = c1 * (vx[[ky, kx]] - vx[[ky, kx - 1]]) - 
							c2 * (vx[[ky, kx + 1]] - vx[[ky, kx - 2]]);
					
						vyy = c3 * (vy[[ky, kx]] - vy[[ky - 1, kx]]) - 
							c4 * (vy[[ky + 1, kx]] - vy[[ky - 2, kx]]);

						(*Update velocity*)
						p[[ky, kx]] = p[[ky, kx]] - l[[ky, kx]] * dt * (vxx + vyy);,

						{ky, 5, ny - 5, 1}
					];,
					{kx, 5, nx - 5, 1}
				]; 
				
				(*Return*)
				p, 
				
				{q, signal}
			]
		], 
		
		CompilationTarget -> target, 
		CompilationOptions -> {"ExpressionOptimization" -> False}, 
		RuntimeOptions -> "Speed"
	];


(* ::Section:: *)
(*LWaveSimulation*)


LWaveSimulation::notimpl = 
"Method `1` not implemented"


Options[LWaveSimulation] = 
	{"Method" -> "Compiled", "CompilationTarget" -> "MVM"}


LWaveSimulation[lv: {{__Real}..}, rho: {{__Real}..}, signal: {__Real}, {xs_Integer, ys_Integer}, OptionsPattern[]] /; 
	0 < xs <= Length[lv[[1]]] && 0 < ys <= Length[lv] && Dimensions[lv] === Dimensions[rho] && MatrixQ[lv] := 
	Switch[OptionValue["Method"], 
		"Compiled", 
		LWaveSimulationCompiled[OptionValue["CompilationTarget"]][lv, rho, signal, xs, ys], 
		
		_, 
		Message[LWaveSimulation::notimpl, OptionValue["Method"]]; Null
	]


(* ::Section:: *)
(*End private context*)


End[]


(* ::Section:: *)
(*Protection*)


Protect["`*"]


(* ::Section:: *)
(*End package*)


EndPackage[] (*GPNTools`FDWaveSim`*)
