(* :Title: BlockUnprotect *)
(* :Context: BlockUnprotect` *)
(* :Author: Mark A. Caprio, Center for Theoretical Physics, 
  Yale University *)
(* :Summary: 
      Dynamic scoping structure for options associated with a symbol. *)
(* :Copyright: Copyright 2005, Mark A. Caprio *)
(* :Package Version: 0.0 *)
(* :Mathematica Version: 4.0 *)
(* :History:
      Adapted from BlockOptions December 1, 2005.
  *)



BeginPackage["BlockUnprotect`"];

Unprotect[Evaluate[$Context<>"*"]];





BlockUnprotect::usage="BlockUnprotect[{symbol1,...},body] temporarily unprotects the given symbols while body is evaluated and then restores their prior protection status.";



BlockUnprotect::needlist="The first argument of BlockUnprotect must be a list of symbols.";
BlockUnprotect::numargs="BlockUnprotect must be called with exactly two arguments.";



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





SetAttributes[BlockUnprotect,HoldRest];

BlockUnprotect[x_/;!MatchQ[x,{___Symbol}],_]:=
    Message[BlockUnprotect::needlist];
BlockUnprotect[_]:=Message[BlockUnprotect::numargs];
BlockUnprotect[_,_,__]:=Message[BlockUnprotect::numargs];

BlockUnprotect[IdentifierList:{___Symbol},Body_]:=Module[
      {Identifier,EvaluatedBody,IsProtected,Aborted},
      AbortProtect[
        
        (* Clear protection *)
        DoForEach[
          
          IsProtected[Identifier]=
            MemberQ[Attributes[Evaluate[Identifier]],Protected];
          If[IsProtected[Identifier],Unprotect[Evaluate[Identifier]]],
          {Identifier,IdentifierList}
          ];
        
        (* Evaluate body *)
        Aborted=False;
        CheckAbort[
          EvaluatedBody=Body,
          Aborted=True
          ];
        
        (* Restore protection *)
        DoForEach[
          If[IsProtected[Identifier],Protect[Evaluate[Identifier]]],
          {Identifier,IdentifierList}
          ];
        ];
      
      (* Return value *)
      (* Passes through abort, 
        and also explicitly returns $Aborted in case Abort[] is suppressed *)
      
      If[Aborted,Abort[];$Aborted,EvaluatedBody]
      ];





End[];



Protect[Evaluate[$Context<>"*"]];
EndPackage[];





































