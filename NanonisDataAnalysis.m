(* Nanonis Data Analysis Package *)

BeginPackage[ "NanonisDataAnalysis`"]

Print["Let's begin!"]

NanonisOpener::usage = 
	"NanonisOpener identifies the type of nanonis file from its extension, reads it and outputs a Matematica association with the file properties" 

NanonisFileParser::usage = 
	"NanonisFileParser[filepath] identifies and reads the file, outputs a simple association: including only the Header and Data keys"

PixelData::usage = "PixelData[dataset,pixelVariableName]Extracts the data from a 3ds file and parses it into data per pixel, stored as an association in the user defined variable name"

MakeMap::usage = "MakeMap[dataset,mapVariableName]Extracts the data from a 3ds file and parses it into data per pixel, stored as an association in the user defined variable name"

MapFilter::usage = "MapFilter[map, varStored, filterRadius] filters out a map using a filtering function (Gaussian is default)"

MakeCoord::usage = "MakeCoord[x,y] quickly makes coordinate pairs from two lists."

TakeMap::usage = "Extracts the map from the association with the key closed to the selected value."

MapPlot::usage = "Plots the parsed maps"

ScanParser::usage = "Takes in scan file and parses the data in it into its different channels"

ScanViewer::usage = "ScanViewer[scanFile] allows you to view Scan Files"

ScanTo3DCoordinates::usage = "Take a scan file and transofm it into a list of 3D points"

PlaneSubstract::usage ="Substract a Plane"

TensorToPixel::usage="TensorToPixel[array,pixVar,dataset] takes a tensor array and transforms it into a set of pixel points"

FourierTransformShift::usage = "Shifts a Fourier Matrix"

InverseFourierTransformShift::usage = "Inverse operation to InverseFourierTransformShift"

BresenhamLine::usage = "BresenhamLine[pt0,pt1] Gives the coordinate between the integer points pt0 and pt1"

Begin["Private`"]

MakeCoord[x_,y_]:=Transpose[{x,y}];

NanonisFileParser[filepath_]:= Module[{fileExtension, fileBaseName,
	workingStream, headerEndPosition, data, header, dataOutput
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
		Find[workingStream, ":SCANIT_END:"];(*Finds the position the the end of the header binarywise*)
		headerEndPosition = StreamPosition[workingStream];(*Passes the bit position of the header delimiter*)
		SetStreamPosition[workingStream, 
		 headerEndPosition +  5];(*Sets stream position after the stop bits (2) where the data \
		file begins*)
		data = 
		 BinaryReadList[workingStream, "Real32", 
		  ByteOrdering -> 
		   1];(*Reads the rest of the file after the divinding bits and \
		header. The data is read as a 32 byte floating point number with a \
		big endian ordering*)
		SetStreamPosition[workingStream, 0];(*Sets the stream position back \
		to the beginning*)
		header = StringJoin[
		  ReadList[workingStream, "Character", 
		   headerEndPosition]];(*Reads the header as a string and stores it into \
		the header variable*)
		,
		If[fileExtension == "dat",
		Print["Importing point spectroscopy file..."];
		,
		Print["File not supported"]
				];
			];
		];
	Close[workingStream];(*Closes the file stream*)
	dataOutput=Association["Header"-> header,"Data" -> data];
	Print["...Done"];
	dataOutput
]

HeaderParser3DS[header_] := 
 Module[{fullHeader, listHeader, isContainsVar, 
   headerLineVar, parsedHeader,headerVars},
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

HeaderParserSXM[header_] := 
 Module[{splitHeader, headerLength, scanHeaderKeys, scanHeaderValues,parsedHeader},
  splitHeader = 
   StringReplace[
    StringSplit[header, 
     RegularExpression["(:\\w+[[:ascii:]]*?:)"] -> "$0"], "\n" -> ""];
  headerLength = Length[splitHeader];
  scanHeaderKeys = 
    DeleteCases[
     StringReplace[Take[splitHeader, {1, headerLength, 2}], 
      ":" -> ""], "SCANIT_END"];
  scanHeaderValues = Take[splitHeader, {2, Length[splitHeader], 2}];
  parsedHeader = AssociationThread[scanHeaderKeys, scanHeaderValues];
  parsedHeader
  ]

NanonisOpener[filepath_] := Module[{datafile, data, parsedHeader,fileExtension},
  datafile = NanonisFileParser[filepath];
  fileExtension = FileExtension[filepath];
  If[fileExtension=="3ds",
  	parsedHeader = HeaderParser3DS[datafile["Header"]],
  	If[fileExtension=="sxm",
  	parsedHeader=HeaderParserSXM[datafile["Header"]]
  		]
  	];
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
      Table[dIdVTensor[[i, j]][[k]], {i, 1, pixNy}, {j, pixNx}];
    ];
   mapVarName = 
    AssociationThread[bias, Table[maps[k], {k, 1, sweepPoints}]];
   ,
   Print["No LockinX data or the channel name is different."]
   ];
  ];

TakeMap[value_, map_] := map[First[Nearest[Keys[map], value]]];

MapPlot[map_, OptionsPattern[]] := DynamicModule[{minVal,maxVal,gradientColors},
   minVal = Min[Values[map]]; maxVal = Max[Values[map]];
   gradientColors = 
    Thread[ColorData["Gradients"] -> ColorData["Gradients", "Image"]];
  Manipulate[
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

scanReshape[data_,pixNx_,pixNy_] := ArrayReshape[data, {pixNy, pixNx}];
fwdDataArray[data_,isUp_,pixNx_,pixNy_] := If[isUp, Reverse[scanReshape[data,pixNx,pixNy]],
   scanReshape[data,pixNx,pixNy]
   ];
bwdDataArray[data_,isUp_,pixNx_,pixNy_] := If[isUp, Reverse[Reverse[scanReshape[data,pixNx,pixNy]], 2],
   Reverse[scanReshape[data,pixNx,pixNy], 2]
   ];

ScanParser[scanFile_] := 
 Module[{scanHeader, pixNx, pixNy, pixRx, pixRy, pixTotal, 
   scanChannels, scanChannelsN, entriesPerChannel, 
   scanDataChannelPartitioned, isUp, scanData, i, channelName, 
   channelData, channelFwd, channelBwd,rangeX,rangeY,scanCoordMesh},
  scanHeader = scanFile["Header"];
  {pixNx, pixNy} = 
   ToExpression[
    StringCases[scanHeader["SCAN_PIXELS"], RegularExpression["\\d+"]]];
  {pixRx, pixRy} = 
   ToExpression[
     StringReplace[StringSplit[scanHeader["SCAN_RANGE"], "      "], 
      "E" :> "*^"]]*10^9;
  pixTotal = Times @@ {pixNx, pixNy};
  rangeX = Subdivide[0, pixRx, pixNx - 1];
  rangeY = Subdivide[0, pixRy, pixNy - 1];
  scanCoordMesh = 
  Table[{(rangeX)[[j]], (Reverse@rangeY)[[i]]}, {i, 1, pixNx}, {j, 1, 
    pixNy}];
  isUp = StringMatchQ[scanHeader["SCAN_DIR"], "up"];
  scanChannels = 
   StringDelete[
     StringReplace[StringSplit[scanHeader["Scan>channels"], ";"], 
      "(" ~~ _ ~~ ")" -> ""], Whitespace] // Reverse;
  scanChannelsN = Length[scanChannels];
  entriesPerChannel = 
   Quotient[FromDigits@Dimensions[scanFile["Data"]], 
    pixNx*pixNy*scanChannelsN];
  scanDataChannelPartitioned = 
   Partition[scanFile["Data"], {pixNx*pixNy*entriesPerChannel}];
  scanData = 
   Association["Pixels" -> {pixNx, pixNy}, "TotalPixels" -> pixTotal, 
    "Range" -> {pixRx, pixRy},"Range Lists"->{rangeX,rangeY},"Coordinate Mesh"->scanCoordMesh,"Channels"->scanChannels];
  For[i = 1, i <= Length[scanChannels], i++,
   channelName = scanChannels[[i]];
   channelData = 
    Partition[scanDataChannelPartitioned[[i]], {pixNx*pixNy}];
   channelFwd = fwdDataArray[channelData[[1]],isUp,pixNx,pixNy];
   channelBwd = bwdDataArray[channelData[[2]],isUp,pixNx,pixNy];
   AssociateTo[
    scanData, {channelName <> " Fwd" -> channelFwd, 
     channelName <> " Bwd" -> channelBwd}]
   ];
  scanData
  ]

ScanViewer[scanFile_] := DynamicModule[{scanDataset,channel,direction},
  scanDataset = ScanParser[scanFile];
  Manipulate[
   Module[{channels, ticksX, ticksY, pixNx, pixNy, pixRx, pixRy},
    {pixNx, pixNy} = scanDataset["Pixels"];
    {pixRx, pixRy} = scanDataset["Range"];
    ticksX[min_, max_] := 
     Table[{i, Round[(i) pixRx/(pixNx), 0.01]}, {i, min, max, pixRx}];
    ticksY[min_, max_] := 
     Table[{i, Round[(i) pixRy/(pixNy), 0.01]}, {i, min, max, pixRy}];
    ArrayPlot[scanDataset[channel <> direction]
     , ColorFunction -> ColorData["DeepSeaColors"]
     , FrameTicks -> {ticksY, ticksX}, ImageSize -> Full
     ]
    ], {{channel, "Z", "Channel"}, scanDataset["Channels"], 
    ControlType -> RadioButton}, {{direction, " Fwd", 
     "Direction"}, {" Fwd" -> "Forward", " Bwd" -> "Backward"}}
   ]
  ]
  
 ScanTo3DCoordinates[scanDataset_, channel_, flat_: False] := 
 Module[{pixNx, pixNy, rangeX, rangeY, scan3DCoordArray},
  {pixNx, pixNy} = scanDataset["Pixels"];
  {rangeX, rangeY} = scanDataset["Range Lists"];
  scan3DCoordArray = 
   Table[{rangeX[[j]], (Reverse@
        rangeY)[[i]], (scanDataset[channel][[i, j]])*10^9}, {i, pixNy}, {j, pixNx}];
  If[flat, Flatten[scan3DCoordArray, 1], scan3DCoordArray]
  ]
  
 PlaneSubstract[scanFile_, channel_] := 
 Module[{scanDataset, pixNx, pixNy, scan3DCoordFlat, fittedPlaneArray},
  scanDataset = ScanParser[scanFile];
  {pixNx, pixNy} = scanDataset["Pixels"];
  scan3DCoordFlat = ScanTo3DCoordinates[scanDataset, channel, True];
  fittedPlaneArray = 
   Partition[
    Fit[scan3DCoordFlat, {1, x, y}, {x, y}, "PredictedResponse"], 
    pixNx];
 (scanDataset[channel])*10^9 - fittedPlaneArray
  ]
  
Options[MapFilter] = {GaussianFiltering -> True};
MapFilter[maps_, varName_, filterSize_, OptionsPattern[]] := 
 Module[{filterFunct},
  If[OptionValue[GaussianFiltering],
   Print["Gaussian Smoothing Map..."];
   filterFunct[map_] := GaussianFilter[map, filterSize];
   ];
  varName = Map[filterFunct, maps];
  Print["...Done"];
  ]

TensorToPixel[array_, varName_, dataset_] := 
 Module[{pixNx, pixNy, i, j},
  {pixNx, pixNy} = 
   ToExpression[StringSplit[dataset["Header"]["Grid dim"], "x"]];
  For[j = 0, j < pixNy, j++,
   For[i = 1, i <= pixNx, i++,
    varName[i + (pixNx*j)] = array[[All, pixNy - j, i]];
    ]
   ]
  ]
  
FourierTransformShift[dat_?ArrayQ, k : (_Integer?Positive | All) : All] :=
         Module[{dims = Dimensions[dat]}, 
                RotateRight[dat, If[k === All, Quotient[dims, 2], 
                                    Quotient[dims[[k]], 2] UnitVector[Length[dims], k]]]]

InverseFourierTransformShift[dat_?ArrayQ, k : (_Integer?Positive | All) : All] := 
          Module[{dims = Dimensions[dat]}, 
                 RotateRight[dat, If[k === All, Ceiling[dims/2], 
                                     Ceiling[dims[[k]]/2] UnitVector[Length[dims], k]]]]

lowLine[pt0_,pt1_]:=Module[{x0,x1,y0,y1,yi,dx,dy,d,y,i,xCoords,yCoords},
{x0,y0}=pt0;{x1,y1}=pt1;
dx=x1-x0;dy=y1-y0;yi=1;
If[dy<0,yi=-1;dy=-dy];
d=2*dy-dx;
y=y0;
xCoords=Table[x,{x,x0,x1,1}];
yCoords={};
For[i=1,i<=Length[xCoords],i++,
yCoords=Append[yCoords,y];
If[d>0,y=y+yi;d=d-2*dx];
d=d+2*dy;
];
Transpose@{xCoords,yCoords}
]

highLine[pt0_,pt1_]:=Module[{x0,x1,y0,y1,yi,dx,dy,d,y,i,x,xi,xCoords,yCoords},
{x0,y0}=pt0;{x1,y1}=pt1;
dx=x1-x0;dy=y1-y0;xi=1;
If[dx<0,xi=-1;dx=-dx];
d=2*dx-dy;
x=x0;
yCoords=Table[y,{y,y0,y1,1}];
xCoords={};
For[i=1,i<=Length[yCoords],i++,
xCoords=Append[xCoords,x];
If[d>0,x=x+xi;d=d-2*dy];
d=d+2*dx;
];
Transpose@{xCoords,yCoords}
]

BresenhamLine[pt0_,pt1_]:=Module[{x0,y0,x1,y1},
{x0,y0}=pt0;{x1,y1}=pt1;
If[Abs[y1-y0]<Abs[x1-x0],
If[x0>x1,lowLine[pt1,pt0],lowLine[pt0,pt1]],
If[y0>y1,highLine[pt1,pt0],highLine[pt0,pt1]]
]
]



End[]

EndPackage[]

