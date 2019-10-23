(* ::Package:: *)

(*Peak Finding Package*)

BeginPackage["dIdVPeakFind`"]

dIdVInterpolation::usage = "dIdVInterpolation[pixBias,pixdIdV,pixTotal,dIdVInterVar,dIdVDInterVar] takes in the bias and spectrum data (pixBias,pixdIdV) at each pixel as well as the total number of pixels and creates an interpolating function of the curve and its derivative and stores them in the user defined variables dIdVInterVar and dIdVDInterVar"
PeakFind::usage = "PeakFind[Bias,pixTotal,dIdVDInterVar,minTrialBias,maxTrialBias,minRootSearch,maxRootSearch,DThresh] Takes in the bias, total number of pixels and the derivative of the dIdV interpolated curve (Bias,pixTotal,dIdVDInterVar) along with the user defined variables of the minimum and maximum bias trial ranges (minTrialBias, maxTrialBias) where one expects the peaks to be in, and the bias range that the root finding algorithm will walk to try to find a peak (minRootSearch, maxRootSearch) as well as the threshold to use on the second derivative of the dIdV curve to find the \[OpenCurlyDoubleQuote]shoulder-like\[CloseCurlyDoubleQuote] features. Stores the output as lists under the variable names roots[pixNum] and roots2D[pixNum] "
CheckPHSymmetry::usage = "CheckPHSymmetry[\[Delta],\[Delta]sh,pixTotal,rootList,root2DList] Takes the energy window radius inside where to look for PH symmetry for the peaks and shoulder like features (\[Delta],\[Delta]sh) in units of \[CapitalDelta]E in the rootList and shoulder lists for each pixel. Outputs the results as a list per pixel writing it to the variable \[CapitalDelta]"

Begin["Private`"]

dIdVInterpolation[pixBias_,pixdIdV_,pixTotal_,dIdVInterVar_,dIdVDInterVar_]:=Module[{i},
For[i=1,i<=pixTotal,i++,
dIdVInterVar[i] =Interpolation[MakeCoord[pixBias[i],pixdIdV[i]]];
dIdVDInterVar[i]=Interpolation[MakeCoord[Drop[pixBias[i],-1],Differences[pixdIdV[i]]/Differences[pixBias[i]]]]
]
]

PeakFind[Bias_,pixTotal_,dIdVDInterVar_,minTrialBias_,maxTrialBias_,minRootSearch_,maxRootSearch_,DThresh_]:=DynamicModule[{n,i,rootP,rootN,root2Dp,root2Dn,rootSearchBiasP,rootSearchBiasN},
CreateDialog[Column[{"Finding Peaks...",Row[{ProgressIndicator[Dynamic[n],{1,pixTotal}]," ",Dynamic[Round@N[100 n/pixTotal]]," %"}]}],WindowSize->Fit];(**)
Quiet[
For[n=1,n<=pixTotal,n++,
roots[n]={};
roots2D[n]={};
rootSearchBiasP=Cases[Bias[n],x_/;x>Abs[minTrialBias]&&x<Abs[maxTrialBias]];
rootSearchBiasN=Cases[Bias[n],x_/;x<-Abs[minTrialBias]&&x>-Abs[maxTrialBias]];
For[i=1,i<=Length[rootSearchBiasP],i++,
rootP=FindRoot[dIdVDInterVar[n][e],{e,rootSearchBiasP[[i]],Abs[minRootSearch],Abs[maxRootSearch]},AccuracyGoal->1][[1,2]];
If[Negative[D[dIdVDInterVar[n][e],e]/.e->rootP],
AppendTo[roots[n],rootP];
];
rootN = FindRoot[dIdVDInterVar[n][e],{e,rootSearchBiasN[[i]],-Abs[maxRootSearch],-Abs[minRootSearch]},AccuracyGoal->1][[1,2]];
If[Negative[D[dIdVDInterVar[n][e],e]/.e->rootN],
AppendTo[roots[n],rootN];
];
root2Dp=FindRoot[D[dIdVDInterVar[n][e],e],{e,rootSearchBiasP[[i]],Abs[minRootSearch],Abs[maxRootSearch]}][[1,2]];
If[Positive[D[dIdVDInterVar[n][e],{e,2}]/.e->root2Dp]&&(0<=dIdVDInterVar[n][root2Dp]<=Abs[DThresh]),
AppendTo[roots2D[n],root2Dp]
];
root2Dn=FindRoot[D[dIdVDInterVar[n][e],e],{e,rootSearchBiasN[[i]],-Abs[maxRootSearch],-Abs[minRootSearch]}][[1,2]];
If[Negative[D[dIdVDInterVar[n][e],{e,2}]/.e->root2Dp]&&(0>=dIdVDInterVar[n][root2Dp]>=-Abs[DThresh]),
AppendTo[roots2D[n],root2Dn]
];
]
]
]
]

CheckPHSymmetry[\[CapitalDelta]E_,\[Delta]_,\[Delta]sh_,pixTotal_,rootList_,root2DList_]:=Module[{n,negRoots,posRoots,shortList,shortListPHEnergies,longList,longListPHEnergies,sortedListbyLength,found\[CapitalDelta],tmp\[CapitalDelta],allroots,neg2DRoots,pos2DRoots,short2DList,long2DList},
For[n=1,n<=pixTotal,n++,
\[CapitalDelta][n]={};
tmp\[CapitalDelta]={};
(*Peak Routine*)
negRoots[n]=Abs[Cases[rootList[n],x_/;Negative[x]]];
posRoots[n]=Abs[Cases[rootList[n],x_/;Positive[x]]];
If[Length@negRoots[n]!=0&&Length@posRoots[n]!=0,
sortedListbyLength = SortBy[{negRoots[n],posRoots[n]},Length];
shortList=sortedListbyLength[[1]];
longList=sortedListbyLength[[-1]];
longListPHEnergies = Flatten[NearestTo[shortList,{1,Round[\[CapitalDelta]E*\[Delta],0.01]}][longList]];
shortListPHEnergies=Flatten[Nearest[shortList,longListPHEnergies,{1,Round[\[CapitalDelta]E*\[Delta],0.01]}]];
found\[CapitalDelta] = DeleteDuplicates[Round[Mean[{longListPHEnergies,shortListPHEnergies}],0.0001]];
AppendTo[tmp\[CapitalDelta],found\[CapitalDelta]];
];
(*Shoulder Routine*)
allroots[n]=Join[rootList[n],root2DList[n]];
neg2DRoots[n]=Abs[Cases[allroots[n],x_/;Negative[x]]];
pos2DRoots[n]=Abs[Cases[allroots[n],x_/;Positive[x]]];
If[Length@neg2DRoots[n]!=0&&Length@pos2DRoots[n]!=0,
sortedListbyLength = SortBy[{neg2DRoots[n],pos2DRoots[n]},Length];
short2DList=sortedListbyLength[[1]];
long2DList=sortedListbyLength[[-1]];
longListPHEnergies = Flatten[NearestTo[short2DList,{1,Round[\[CapitalDelta]E*\[Delta]sh,0.01]}][long2DList]];
shortListPHEnergies=Flatten[Nearest[short2DList,longListPHEnergies,{1,Round[\[CapitalDelta]E*\[Delta]sh,0.01]}]];
found\[CapitalDelta] = DeleteDuplicates[Round[Mean[{longListPHEnergies,shortListPHEnergies}],0.0001]];
AppendTo[tmp\[CapitalDelta],found\[CapitalDelta]];
];
AppendTo[\[CapitalDelta][n],DeleteDuplicates[Flatten[tmp\[CapitalDelta]]]];
\[CapitalDelta][n]=Flatten[\[CapitalDelta][n]];
];
]

End[]

EndPackage[]
