(* ::Package:: *)

(* Nanonis Data Analysis Package *)

BeginPackage[ "NanonisDataAnalysis`"]

Print["Let's begin!"]

NanonisOpener::usage = 
	"NanonisOpener identifies the type of nanonis file from its extension, reads it and outputs a Matematica association with the file properties" 

NanonisFileParser::usage = 
	"NanonisFileParser[filepath] identifies and reads the file, outputs a simple association: including only the Header and Data keys"

PixelData::usage = "PixelData[dataset,pixelVariableName]Extracts the data from a 3ds file and parses it into data per pixel, stored as an association in the user defined variable name"

MakeMap::usage = "MakeMap[dataset,mapVariableName]Extracts the data from a 3ds file and parses it into data per pixel, stored as an association in the user defined variable name"

MakeCoord::usage = "MakeCoord[x,y] quickly makes coordinate pairs from two lists."

TakeMap::usage = "Extracts the map from the association with the key closed to the selected value."

MapPlot::usage = "Plots a Grid Spectroscopy Map"


Begin["Private`"]

MakeCoord[x_,y_]:=Transpose[{x,y}];

NanonisFileParser[filepath_]:= Module[{fileExtension, fileBaseName,
	workingStream, headerEndPosition, data, header
	},
	(*Identify File Name and Format*)
	fileExtension = FileExtension[filepath];
	fileBaseName = FileBaseName[filepath];
	workingStream = OpenRead[filepath, BinaryFormat -> True];(*Opens a binary stream into the selected file*)
	If[fileExtension == "3ds",
		Print["Importing 3ds file..."];
		Find[workingStream, ":HEADER_END:"];(*Finds the position the the end \
		of the header binarywise*)
		headerEndPosition = 
		 StreamPosition[workingStream];(*Passes the bit position of the header delimiter*)
		SetStreamPosition[workingStream, headerEndPosition + 2];(*Sets stream position after the stop bits (2) where the data \
		file begins*)
		data = BinaryReadList[workingStream, "Real32", ByteOrdering -> 1];(*Reads the rest of the file after the divinding bits and \
		header. The data is read as a 32 byte floating point number with a \
		big endian ordering*)
		SetStreamPosition[workingStream, 0];(*Sets the stream position back \
		to the beginning*)
		header = StringJoin[ReadList[workingStream, "Character", headerEndPosition]];(*Reads the header as a string and stores it into \
		the header variable*)
		,
		If[fileExtension == "sxm",
		Print["Importing scan file..."];
		,
		If[fileExtension == "dat",
		Print["Importing point spectroscopy file..."];
		,
		Print["File not supported"]
				];
			];
		];
	Close[workingStream];(*Closes the file stream*)
	Association["Header"-> header,"Data" -> data]
]

HeaderParser[header_] := 
 Module[{fullHeader, listHeader, isContainsVar, 
   headerLineVar, parsedHeader},
  fullHeader = {};
  headerVars = {"Grid dim", "Grid settings", "Sweep Signal", 
    "Fixed parameters", "Experiment parameters", "# Parameters", 
    "Experiment size", "Points", "Channels", "Experiment", 
    "Delay before measuring", "Start time", "End time", "Comment"};
  For[i = 1, i <= Length[headerVars], i++,
   If[StringContainsQ[header, headerVars[[i]]],
     listHeader = StringSplit[header, "\n"];
     isContainsVar = StringContainsQ[listHeader, headerVars[[i]]];
     headerLineVar = 
      StringSplit[
       StringReplace[Pick[listHeader, isContainsVar][[1]], 
        "\"" -> ""], "="];
     AppendTo[fullHeader, headerLineVar]
     ];
   ];
  parsedHeader = 
   Association[
    Thread[fullHeader[[All, 1]] -> 
      StringReplace[fullHeader[[All, 2]], Whitespace -> ""]]];
  parsedHeader
  ]

NanonisOpener[filepath_] := Module[{datafile, data, parsedHeader},
  datafile = NanonisFileParser[filepath];
  parsedHeader = HeaderParser[datafile["Header"]];
  data = datafile["Data"];
  Association[{"Header" -> parsedHeader, "Data" -> data}]
  ]
 
 stringPreprocessing [string_] := 
 StringSplit[StringReplace[string, "(" ~~ _ ~~ ")" -> ""], ";"]
  
 PixelData[dataset_, pixVarName_] := 
 Module[{pixNx, pixNy, pixN, fixedParams, expParams, allParams, 
   paramsN, channels, channelsN, sweepPoints,
   pixArray, i, pixParams, pixPartionedData},
  {pixNx, pixNy} = 
   ToExpression[StringSplit[dataset["Header"]["Grid dim"], "x"]];
  pixN = pixNx*pixNy;
  fixedParams = 
   stringPreprocessing[dataset["Header"]["Fixed parameters"]];
  expParams = 
   stringPreprocessing[dataset["Header"]["Experiment parameters"]];
  allParams = Join[fixedParams, expParams];
  paramsN = Length[allParams];
  channels = stringPreprocessing[dataset["Header"]["Channels"]];
  channelsN = Length[channels];
  sweepPoints = ToExpression[dataset["Header"]["Points"]];
  pixArray = 
   Partition[dataset["Data"], channelsN*sweepPoints + paramsN];
  For[i = 1, i <= pixN, i++,
   pixParams = Take[pixArray[[i]], paramsN];
   pixPartionedData = 
    Partition[Drop[pixArray[[i]], paramsN], sweepPoints];
   pixVarName[i] = AssociationThread[allParams, pixParams];
   AssociateTo[
    pixVarName[i], {AssociationThread[channels, pixPartionedData]}];
   ];
  ]
 
 MakeMap[dataset_, mapVarName_] := 
 Module[{pixData, pixNx, pixNy, sweepPoints, bias, sweepStart, 
   sweepEnd,
   dIdVTensor, maps,i,j,k},
  PixelData[dataset, pixData];
  {pixNx, pixNy} = 
   ToExpression[StringSplit[dataset["Header"]["Grid dim"], "x"]];
  sweepPoints = ToExpression[dataset["Header"]["Points"]];
  If[KeyMemberQ[pixData[1], "Bias"],
   bias = ToExpression[(pixData[1]["Bias"])]*10^3;
   Print[bias];
   ,
   Print["No Bias channel saved, bias will be infered from the \
initial and final sweep setting and number of points"];
   {sweepStart, sweepEnd} = 
    ToExpression[ {pixData[1]["SweepStart"], pixData[1]["SweepEnd"]}];
   bias = Subdivide[sweepStart, sweepEnd, sweepPoints - 1]*10^3;
   ];
  If[KeyMemberQ[pixData[1], "LockinX"],
   dIdVTensor = 
    Table[pixData[i + (pixNx*(pixNy - j))]["LockinX"], {j, 1, 
      pixNy}, {i, 1, pixNx}];
   For[k = 1, k <= sweepPoints, k++,
    maps[k] = 
      Table[dIdVTensor[[i, j]][[k]], {i, 1, pixNx}, {j, pixNy}];
    ];
   mapVarName = 
    AssociationThread[bias, Table[maps[k], {k, 1, sweepPoints}]];
   ,
   Print["No LockinX data or the channel name is different."]
   ];
  ];

TakeMap[value_, map_] := map[First[Nearest[Keys[map], value]]];

MapPlot[map_, OptionsPattern[]] := DynamicModule[{minVal,maxVal},
  Manipulate[
   minVal = Min[Values[map]]; maxVal = Max[Values[map]];
   gradientColors = 
    Thread[ColorData["Gradients"] -> ColorData["Gradients", "Image"]];
   ArrayPlot[TakeMap[Bias, map],
    PlotLabel -> ("Bias: " <> ToString[Bias]),
    ColorFunction -> plotColor,
    PlotLegends -> True,
    If[OptionValue[DynamicScaling] == False, 
     PlotRange -> {minVal, maxVal}, PlotRange -> Automatic]
    ],
   {{Bias,Min[Keys[map]],"Bias"}, Keys[map],Manipulator}, {{plotColor, "DeepSeaColors","Gradient"}, 
    gradientColors,ControlType->PopupMenu},
   AppearanceElements -> All, ContinuousAction->True
   ]
]
Options[MapPlot] = {DynamicScaling -> False};


End[]

EndPackage[]



