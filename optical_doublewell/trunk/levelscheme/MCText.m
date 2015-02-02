(* :Title: MCText *)
(* :Context: MCText` *)
(* :Author: Mark A. Caprio, Center for Theoretical Physics, 
  Yale University *)
(* :Summary: Miscellaneous text formatting utilities. *)
(* :Copyright: Copyright 2005, Mark A. Caprio *)
(* :Package Version: 0.0 *)
(* :Mathematica Version: 4.0 *)
(* :History:
      Package started December 2003, 
  based upon functions written for LevelScheme.
      Distributed with the LevelScheme package.
      Rev. January 25, 2006.
  *)



BeginPackage["MCText`"];

Unprotect[Evaluate[$Context<>"*"]];





textup::usage="textup[text] gives non-italic text.";
textsl::usage="textsl[text] gives slanted text.";
textit::usage="textit[text] gives italic text.";
textmd::usage="textmd[text] gives non-bold text.";
textbf::usage="textbf[text] gives bold text.";
textrm::usage="textrm[text] gives text in the Times font.";
texttt::usage="texttt[text] gives text in the Courier font.";
textsf::usage="textsf[text] gives text in the Helvetica font.";

textsize::usage="textsize[size,text] gives text with the specified point size.";\

textcolor::usage="textcolor[color,text] gives text with the specified color.";\

texttracking::usage="texttracking[tracking,text] gives text with the specified tracking.";\

textfamily::usage="textfamily[family,text] gives text with the specified family.";\

texthidden::usage="texthidden[text] gives invisible text.";

textsubscript::usage="UNDOCUMENTED";
textsuperscript::usage="UNDOCUMENTED";
textsubsuperscript::usage="UNDOCUMENTED";

hspace::usage="hspace[width] produces box of given width in ems, which may be negative.";\


StackText::usage="StackText[alignment,linegap,{line1,...}] produces a multiline label.";\


LabelJP::usage="LabelJP[J,P], or Label[J] for P=+1, produces a level spin label, with rational number fractions converted to typeset fractions.";\

LabelJiP::usage="LabelJiP[J,i,P], or Label[J,i] for P=+1, produces a level spin label, with rational number fractions converted to typeset fractions, with subscript i.";\


NucleusBox::usage="UNDOCUMENTED (Limitation: Left subscript and superscript are left aligned.)";\

NuclearA::usage="UNDOCUMENTED";
NuclearN::usage="UNDOCUMENTED";
NuclearZ::usage="UNDOCUMENTED";

SuperPrimeBox::usage="SuperPrimeBox[x] places x in a SuperscriptBox with superscript prime. (UNDOCUMENTED)";\

SuperPrime::usage=
    "SuperPrime[x] superscripts x with a prime, in analogy to SuperPlus, SuperMinus, etc. -- since uses Superscript, behaves like mathematical expression, including automatic parenthesization (DEPRECATED)";\


MultipletBox::usage="UNDOCUMENTED";
TightRowBox::usage="UNDOCUMENTED";

SolidusFractionBox::usage="SolidusFractionBox[x,y] typesets x/y as a solidus fraction.";\

DiagonalFractionBox::usage="DiagonalFractionBox[x,y] typesets x/y as a diagonal fraction.  The option ColumnSpacings->(spacing) controls the horizontal separation between the elements of the fraction.  The option Baseline->{numeratorheight,slashheight,denominatorheight} controls the vertical positioning of the different elements of the fraction.  KernForSuperscript->(adjustment) introduces a horizontal adjustment to the position of any superscript attached to the fraction.  The default option values produce visually reasonable results for the Times New Roman font.";\

KernForSuperscript::usage="KernForSuperscript is an option for DiagonalFractionBox.";\


Fractionize::usage="Fractionize[expr] converts any expression with nonzero denominator into a FractionBox.";\

SolidusFractionize::usage="SolidusFractionize[expr] converts any expression with nonzero denominator into a SolidusFractionBox.";\

DiagonalFractionize::usage="DiagonalFractionize[expr] converts any expression with nonzero denominator into a DiagonalFractionBox.";\


FractionString::usage="Provides a compact solidus-delimited fraction representation of a rational number, in string form. (UNDOCUMENTED)";\

PiFractionString::usage="Provides a compact solidus-delimited fraction representation of a rational multiple of Pi, in string form. (UNDOCUMENTED)";\


SignString::usage="SignString[x] returns \"+\", \"\", or \"-\", depending upon the value of Sign[x].";\


Sqrtize::usage="Sqrtize[x] formats a fraction which is the square root of a rational number entirely under the radical (e.g., the usual format for Clebsch-Gordan coefficients).";\

Radicalize::usage="Radicalize[x,n] formats a fraction involving nth roots entirely under the radical.";\


SubmatrixEllipsis::usage="SubmatrixEllipsis[m,{rows,cols}] prints a submatrix of m with ellipses appended.";\


AlignmentBox::usage= \
"UNDOCUMENTED  (AlignmentBox options: Any option for GridBox may be used. ColumnWidths is crucial option for alignment across labels. Options[AlignmentBox] overrides GridBox defaults of ColumnAlignments and ColumnSpacings. Any option for StyleBox may be used. Background is the most likely to be needed.)";\

Align::usage="Align is an option for AlignmentBox. (UNDOCUMENTED)";
LaTeXTableEntryValue::usage="LaTeXTableEntryValue[Str] takes a typical LaTeX table entry and converts it to a number, by stripping it of any leading nonnumeric characters, any \"&\", and any trailing error estimate following \"(\" or \"\[PlusMinus]\".   (UNDOCUMENTED)";\


MakeColumnBox::usage="UNDOCUMENTED";
MakeHelpBox::usage="UNDOCUMENTED  (Notebook style sheet should be HelpBrowser.)";\

MakeCaptionBox::usage="UNDOCUMENTED";

MakeButtonForURL::usage="UNDOCUMENTED";

NPrint::usage="NPrint[args] applies N to all numeric arguments before passing the resulting modified arguments to Print. (UNDOCUMENTED)";\


PageBreak::usage="PageBreak[] prints a cell with attribute PageBreakBelow set to True.";



LabelJP::invalidparity="Parity must be +1, -1, or None.";
LabelJiP::invalidparity="Parity must be +1, -1, or None.";
MakeColumnBox::notmult="Number of contents items must be multiple of number of columns.";



Begin["`Private`"];

Needs["Utilities`FilterOptions`"]



textup[x_]:=StyleForm[x,FontSlant\[Rule]"Plain"];
textsl[x_]:=StyleForm[x,FontSlant\[Rule]"Oblique"];
textit[x_]:=StyleForm[x,FontSlant\[Rule]"Italic"];
textmd[x_]:=StyleForm[x,FontWeight\[Rule]"Plain"];
textbf[x_]:=StyleForm[x,FontWeight\[Rule]"Bold"];
textrm[x_]:=StyleForm[x,FontFamily\[Rule]"Times"];
texttt[x_]:=StyleForm[x,FontFamily\[Rule]"Courier"];
textsf[x_]:=StyleForm[x,FontFamily\[Rule]"Helvetica"];

hspace[Lems_]:=AdjustmentBox["",BoxMargins\[Rule]{{0,Lems},{0,0}}];

textsize[s_,x_]:=StyleForm[x,FontSize\[Rule]s];
textcolor[c_,x_]:=StyleForm[x,FontColor\[Rule]c];
texttracking[t_,x_]:=StyleForm[x,FontTracking\[Rule]t];
textfamily[f_,x_]:=StyleForm[x,FontFamily\[Rule]f];
texthidden[x_]:=StyleForm[x,ShowContents\[Rule]False];

textsubscript[x_]:=SubscriptBox["",x];
textsuperscript[y_]:=SuperscriptBox["",y];
textsubsuperscript[x_,y_]:=SubsuperscriptBox["",x,y];



Options[textit]={hspace\[Rule]0};
textit[x_,Opts___]:=Module[
      {FullOpts=Flatten[{Opts,Options[textit]}]},
      AdjustmentBox[
        StyleForm[x,FontSlant\[Rule]"Italic"],
        BoxMargins\[Rule]{{0,hspace/.FullOpts},{0,0}}
        ]
      ];



StackText[Alignment_,Spacing_,Lines_List,Opts___?OptionQ]:=
    GridBox[{#}&/@Lines,ColumnAlignments\[Rule]Alignment,
      RowSpacings\[Rule]Spacing,Opts];

SuperPrimeBox[x_]:=SuperscriptBox[x,"\[Prime]"];
SuperPrime[x_]:=Superscript[x,"\[Prime]"];

Options[NucleusBox]={NuclearA->"",NuclearZ->"",NuclearN->""};
NucleusBox[Element_,Opts___?OptionQ]:=Module[
      {FullOpts=Flatten[{Opts,Options[NucleusBox]}]},
      RowBox[{
          SubsuperscriptBox["",NuclearZ/.FullOpts,NuclearA/.FullOpts],
          hspace[-0.2],
          SubsuperscriptBox[Element,NuclearN/.FullOpts,""]  (* 
            to match subsuperscript on left for alignment *)
          }]
      ];





LabelJiP[J_,i_,P_:+1]:=SubsuperscriptBox[
      J/.{Rational\[Rule]FractionBox},
      i,
      Switch[P,+1,"+",-1,"-",None,"",_,Message[LabelJP::invalidparity];""]
      ];
LabelJP[J_,P_:+1]:=SuperscriptBox[
      J/.{Rational\[Rule]FractionBox},
      Switch[P,+1,"+",-1,"-",None,"",_,Message[LabelJP::invalidparity];""]
      ];









MultipletBox[ValueSeq__]:=Module[
      {ValueList={ValueSeq}},
      RowBox[Join[
          {"("},
          {First[ValueList]},
          Flatten[({",",hspace[-0.2],#})&/@Rest[ValueList],1] ,
          {")"}
          ]]
      ];





TightRowBox[Row_List]:=GridBox[{Row},ColumnSpacings\[Rule]0];















SolidusFractionBox[x_,y_,Opts___?OptionQ]:=
    GridBox[{{DisplayForm[x],"/",DisplayForm[y]}},ColumnSpacings\[Rule]-0.1,
      GridBaseline\[Rule]{Baseline,{1,3}}];



Options[DiagonalFractionBox]={ColumnSpacings\[Rule]-0.1,
      Baseline\[Rule]{0.5,0.3,0.0},KernForSuperscript\[Rule]-0.15};
DiagonalFractionBox[x_,y_,Opts___?OptionQ]:=Module[
      {FullOpts=Flatten[{Opts,Options[DiagonalFractionBox]}]},
      TagBox[
        StyleBox[
          GridBox[
            {{
                
                AdjustmentBox[SubscriptBox["",DisplayForm[x]],
                  BoxBaselineShift\[Rule]-((Baseline/.FullOpts)[[1]])],
                
                AdjustmentBox[
                  SubscriptBox["",StyleForm["/",FontSlant->"Oblique"]],
                  BoxBaselineShift\[Rule]-((Baseline/.FullOpts)[[2]])],
                
                AdjustmentBox[SubscriptBox["",DisplayForm[y]],
                  BoxBaselineShift\[Rule]-((Baseline/.FullOpts)[[3]])]
                }},
            ColumnSpacings\[Rule](ColumnSpacings/.FullOpts),
            GridBaseline\[Rule]{Bottom,{1,3}}
            ],
          ScriptBaselineShifts\[Rule]{0,0}
          ],
        DiagonalFractionBox[(KernForSuperscript/.FullOpts)]
        ]
      ];



Unprotect[TagBox];
TagBox/:SuperscriptBox[x:TagBox[_,DiagonalFractionBox[Adjustment_]],n_]:=
    SuperscriptBox[AdjustmentBox[x,BoxMargins\[Rule]{{0,Adjustment},{0,0}}],
      n];
Protect[TagBox];



Fractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=
    DisplayForm[x];
Fractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=
    FractionBox[Numerator[x],Denominator[x],Opts];
f:Fractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];

SolidusFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=
    DisplayForm[x];
SolidusFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=
    SolidusFractionBox[Numerator[x],Denominator[x],Opts];
f:SolidusFractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];

DiagonalFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]==1):=
    DisplayForm[x];
DiagonalFractionize[x_?NumericQ,Opts___?OptionQ]/;(Denominator[x]!=1):=
    DiagonalFractionBox[Numerator[x],Denominator[x],Opts];
f:DiagonalFractionize[x_List,Opts___?OptionQ]:=Thread[Unevaluated[f],List,1];





































































FractionString[x_?NumericQ]:=Module[
      {f,NumeratorString,DenominatorString},
      f=Rationalize[x];
      NumeratorString=ToString[Numerator[f]];
      DenominatorString=ToString[Denominator[f]];
      Which[
        f\[Equal]0,"0",
        Denominator[f]\[Equal]1,StringJoin[NumeratorString],
        Denominator[f]!=1,StringJoin[NumeratorString,"/",DenominatorString]
        ]
      ];



PiFractionString[x_?NumericQ]:=Module[
      {f,NumeratorString,DenominatorString},
      f=Rationalize[x/Pi];
      NumeratorString=If[
          Numerator[f]\[Equal]1,
          "",
          ToString[Numerator[f]]
          ];
      DenominatorString=ToString[Denominator[f]];
      Which[
        f\[Equal]0,"0",
        Denominator[f]\[Equal]1,StringJoin[NumeratorString,"\[Pi]"],
        Denominator[f]!=1,
        StringJoin[NumeratorString,"\[Pi]","/",DenominatorString]
        ]
      ];

















SignString[x_?NumericQ]:=Switch[
      Sign[x],
      +1,"+",
      0,"",
      -1,"-"
      ];



SetAttributes[Sqrtize,Listable];
Sqrtize[x_?NumericQ]:=Which[
      (*RationalQ[x],x,*)
      RationalQ[x^2],Sign[x]*SqrtBox[x^2],
      True,x
      ];
SetAttributes[Radicalize,Listable];
Radicalize[x_?NumericQ,n_Integer]:=Which[
      (*RationalQ[x],x,*)
      RationalQ[x^n],Sign[x]*RadicalBox[Abs[x^n],n],
      True,x
      ];

RationalQ[x_]:=(IntegerQ[x]||(Head[x]===Rational));





























Options[AlignmentBox]={
      AlignmentMarker->"&",
      Align\[Rule]True,
      ColumnAlignments\[Rule]{Right,Left},ColumnSpacings\[Rule]0
      };

BreakString[Separator_,Str_]:=Module[
      {PosnList},
      PosnList=
        Join[{{Null,0}},
          StringPosition[Str,Separator],{{StringLength[Str]+1,Null}}];
      
      Table[StringTake[Str,{PosnList[[i]][[2]]+1,PosnList[[i+1]][[1]]-1}],
        {i,1,Length[PosnList]-1}
        ]
      ];

AlignmentBox[Str_,Opts___?OptionQ]:=Module[
      {
        FullOpts=Flatten[{Opts,Options[AlignmentBox]}]
        },
      
      CheckOption[Align,True|False,FullOpts];
      CheckOption[AlignmentMarker,_String,FullOpts];
      StyleBox[
        If[
          Align/.FullOpts,
          GridBox[
            {BreakString[(AlignmentMarker/.FullOpts),Str]},
            FilterOptions[GridBox,FullOpts]
            ],
          StringReplace[Str,(AlignmentMarker/.FullOpts)->""]
          ],
        FilterOptions[StyleBox,FullOpts]
        ]
      ];









LaTeXTableEntryValue[Value_?NumericQ]:=Value;
LaTeXTableEntryValue[Str_String]:=Module[
      {
        FirstNumericPosn,
        ErrorBarsPosn,
        Value
        },
      
      FirstNumericPosn=
        StringPosition[
            Str,{"+","-","0","1","2","3","4","5","6","7","8","9","."},1][[1,
            1]];
      ErrorBarsPosn=If[
          Length[StringPosition[Str,"(",1]]\[GreaterEqual]1,
          StringPosition[Str,{"(","\[PlusMinus]"},1][[1,1]],
          StringLength[Str]+1
          ];
      Value=ToExpression[
          StringReplace[
            StringTake[Str,{FirstNumericPosn,ErrorBarsPosn-1}],
            {"&"->""}
            ]
          ];
      If[!NumericQ[Value],
        Message[LaTeXTableEntryValue::notnumeric,Str,Value]
        ];
      Value
      ];





SubmatrixEllipsis[m_?MatrixQ,{Rows_,Cols_}]:=
    Append[Transpose[
        Append[Transpose[Take[m,Rows-1,Cols-1]],
          Table["\[CenterEllipsis]",{Cols-1}]]],
      Table["\[CenterEllipsis]",{Rows}]];







MakeButtonForURL[Contents_,URLString_String?(StringMatchQ[#,"http://*"]&),
      Opts___?OptionQ]:=DisplayForm[ButtonBox[
        Contents,
        Opts,
        ButtonData\[RuleDelayed]{URL[URLString],None},ButtonNote->URLString,
        ButtonFunction\[Rule](FrontEndExecute[{FrontEnd`NotebookLocate[#2]}]&)\
,Active\[Rule]True
        ]
      ];





MakeColumnBox[StyleName_String,Columns_Integer,GridBoxOptionsList_List,
    Contents___List]:=Module[
    {GridData},
    
    If[!Mod[Length[{Contents}],Columns\[Equal]0],
      Message[MakeColumnBox::notmult]
      ];
    
    (* Convert each list into a RowBox *)
    GridData=RowBox/@{Contents};
    (* Group entries into rows *)
    GridData=Partition[GridData,Columns];
    (* Convert text styles (created by the LaTeX-
            like directives) into actual boxes *)
    GridData=GridData/.{StyleForm\[Rule]StyleBox};
    
    (* Make cell containing framed grid of table entries *)
    CellPrint[
      Cell[BoxData[
          FormBox[StyleBox[
              FrameBox[
                GridBox[GridData]
                ],
              StyleName,GridBoxOptions\[Rule]GridBoxOptionsList,
              FontFamily->"Times"
              ],DisplayForm]
          ],"Text"]
      ]
    
    ]

MakeCaptionBox[Contents_List]:=Module[
      {ContentsData},
      
      If[!Mod[Length[{Contents}],Columns\[Equal]0],
        Message[MakeColumnBox::notmult]
        ];
      
      (* Convert each list into text data *)
      ContentsData=RowBox[Contents];
      (* Convert text styles (created by the LaTeX-
              like directives) into actual boxes *)
      ContentsData=ContentsData/.{StyleForm\[Rule]StyleBox};
      
      (* Make cell *)
      CellPrint[
        Cell[TextData[ContentsData],"Caption"]
        ]
      
      ];

(* original 2ColumnBox: ColumnWidths\[Rule]{0.31,0.67}*)
MakeHelpBox["SymbolSummary",Contents___List]:=
    MakeColumnBox["2ColumnBox",
      2,{ColumnAlignments\[Rule]{Right,Left},RowLines\[Rule]{False},
        RowSpacings\[Rule]1.5,ColumnWidths\[Rule]{0.31,0.65}},Contents];
MakeHelpBox["SymbolSummaryWideDescription",Contents___List]:=
    MakeColumnBox["2ColumnBox",
      2,{ColumnAlignments\[Rule]{Right,Left},RowLines\[Rule]{False},
        RowSpacings\[Rule]1.5,ColumnWidths\[Rule]{0.26,0.70}},Contents];
(* original 3ColumnBox: {ColumnWidths\[Rule]0.32} *)
MakeHelpBox["OptionSummary",Contents___List]:=
    MakeColumnBox["3ColumnBox",
      3,{ColumnAlignments\[Rule]{Left},RowLines\[Rule]{True,False},
        RowSpacings\[Rule]1.5,ColumnWidths\[Rule]{0.32}},
      {StyleBox["option name","SO10"]},{StyleBox["default value",
          "SO10"]},{""},
      Contents];
MakeHelpBox["OptionSummaryWideDescription",Contents___List]:=
    MakeColumnBox["3ColumnBox",
      3,{ColumnAlignments\[Rule]{Left},RowLines\[Rule]{True,False},
        RowSpacings\[Rule]1.5,ColumnWidths\[Rule]{0.23,0.23,0.5}},
      {StyleBox["option name","SO10"]},{StyleBox["default value",
          "SO10"]},{""},
      Contents];
MakeHelpBox["OptionSummaryWideDescriptionContinuation",Contents___List]:=
    MakeColumnBox["3ColumnBox",
      3,{ColumnAlignments\[Rule]{Left},RowLines\[Rule]{False},
        RowSpacings\[Rule]1.5,ColumnWidths\[Rule]{0.23,0.23,0.5}},
      Contents];











NPrint[ArgSeq___]:=Module[
      {NArgList},
      NArgList=If[NumericQ[#],N[#],#]&/@{ArgSeq};
      Print[Sequence@@NArgList]
      ];







PageBreak:=CellPrint[Cell["",PageBreakBelow\[Rule]True]];



End[];



Protect[Evaluate[$Context<>"*"]];
Unprotect[Evaluate[$Context<>"$*"]];
EndPackage[];





















