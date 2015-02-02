(* :Title: BlockOptions *)
(* :Context: BlockOptions` *)
(* :Author: Mark A. Caprio, Center for Theoretical Physics, 
  Yale University *)
(* :Summary: 
      Dynamic scoping structure for options associated with a symbol. *)
(* :Copyright: Copyright 2005, Mark A. Caprio *)
(* :Package Version: 1.1 *)
(* :Mathematica Version: 4.0 *)
(* :History:
      V1.0, January 22, 2005. MathSource No. 5549.
     V1.1, March 8, 2005.  Attribute HoldAll relaxed to HoldRest.
      V1.2, January 9, 2006.  Added WithOptions.
  *)



BeginPackage["BlockOptions`"];

Unprotect[Evaluate[$Context<>"*"]];





BlockOptions::usage="BlockOptions[{symbol1,...},body] evaluates body with dynamic scoping for Options[symbol1], ..., making local any changes to these options.";\

WithOptions::usage="WithOptions[symbol,{option1->value1,...},body] evaluates body with dynamic scoping for Options[symbol] and with the specified default values for option1, ....";



BlockOptions::needlist="The first argument of BlockOptions must be a list of symbols.";
BlockOptions::numargs="BlockOptions must be called with exactly two arguments.";\

WithOptions::args="WithOptions called with unexpected arguments.";



Begin["`Private`"];



(* Private copy of DoForEach, from package ForEach, version 1.0. *)

SetAttributes[DoForEach,HoldAll];
DoForEach[Expr_,{Var_Symbol,ValueSet_}]:=Module[
      {i},
      
      If[
        Head[ValueSet]=!=List,
        Message[DoForEach::notlist,ValueSet]
        ];
      Do[
        Block[
          {Var=ValueSet[[i]]},
          Expr
          ],
        {i,1,Length[ValueSet]}
        ]
      
      ];
DoForEach[Expr_,{Var_Symbol,CountVar_Symbol,ValueSet_}]:=Module[
      {i},
      
      If[
        Head[ValueSet]=!=List,
        Message[DoForEach::notlist,ValueSet]
        ];
      
      Block[
        {CountVar},
        Do[
          CountVar=i;
          Block[
            {Var=ValueSet[[i]]},
            Expr
            ],
          {i,1,Length[ValueSet]}
          ]
        ]
      
      ];
DoForEach[Expr_,FirstIterator_List,RestIteratorSeq__List]:=
    DoForEach[DoForEach[Expr,RestIteratorSeq],FirstIterator];





SetAttributes[BlockOptions,HoldRest];

BlockOptions[x_/;!MatchQ[x,{___Symbol}],_]:=Message[BlockOptions::needlist];
BlockOptions[_]:=Message[BlockOptions::numargs];
BlockOptions[_,_,__]:=Message[BlockOptions::numargs];

BlockOptions[IdentifierList:{___Symbol},Body_]:=Module[
      {SavedOptions,Identifier,EvaluatedBody,IsProtected,Aborted},
      AbortProtect[
        
        (* Save options *)
        DoForEach[
          SavedOptions[Identifier]=Options[Identifier],
          {Identifier,IdentifierList}
          ];
        
        (* Evaluate body *)
        Aborted=False;
        CheckAbort[
          EvaluatedBody=Body,
          Aborted=True
          ];
        
        (* Restore options *)
        DoForEach[
          IsProtected=MemberQ[Attributes[Evaluate[Identifier]],Protected];
          If[IsProtected,Unprotect[Evaluate[Identifier]]];
          Options[Identifier]=SavedOptions[Identifier];
          If[IsProtected,Protect[Evaluate[Identifier]]],
          {Identifier,IdentifierList}
          ];
        ];
      
      (* Return value *)
      (* Passes through abort, 
        and also explicitly returns $Aborted in case Abort[] is suppressed *)
      
      If[Aborted,Abort[];$Aborted,EvaluatedBody]
      ];



SetAttributes[WithOptions,HoldRest];

WithOptions[x_Symbol,OptList_,Body_]/;MatchQ[OptList,_List?OptionQ]:=
    BlockOptions[
      {x},
      If[Length[OptList]\[GreaterEqual]1,
        SetOptions[x,Sequence@@Flatten[OptList]]];
      Body
      ];
WithOptions[___]/;(Message[WithOptions::args];False):=Null;



End[];



Protect[Evaluate[$Context<>"*"]];
EndPackage[];















































































































