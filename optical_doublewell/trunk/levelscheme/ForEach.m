(* :Title: ForEach *)
(* :Context: ForEach` *)
(* :Author: Mark A. Caprio, Center for Theoretical Physics, 
  Yale University *)
(* :Summary: Control structures for iteration over a set of values. *)
(* :Copyright: Copyright 2005, Mark A. Caprio *)
(* :Package Version: 1.2 *)
(* :Mathematica Version: 4.0 *)
(* :MathSource Number: 5515 *)
(* :History:
      V1.0, January 12, 2005.
      V1.1, January 1, 2006.  Allow multiplet as iteration variable.
      V1.2, September 18, 2006.  SumForEach, ProductForEach.
  *)



BeginPackage["ForEach`"];

Unprotect[Evaluate[$Context<>"*"]];





DoForEach::usage="DoForEach[Expr,{Var,ValueSet},...] or DoForEach[Expr,{Var,CountVar,ValueSet},...] evaluates Expr for each value of Var in the list ValueSet, optionally providing the index (1-based) for each iteration in CountVar.  Last index varies most rapidly.";\

TableForEach::usage="TableForEach[Expr,{Var,ValueSet},...] or TableForEach[Expr,{Var,CountVar,ValueSet},...]  returns a table of Expr, evaluated for each value of Var in the list ValueSet, optionally providing the index (1-based) for each iteration in CountVar.  Last index varies most rapidly.  Lists are nested unless option Flatten is set to true.";\

SumForEach::usage="SumForEach[Expr,{Var,ValueSet},...] or SumForEach[Expr,{Var,CountVar,ValueSet},...]  returns a sum of Expr, evaluated for each value of Var in the list ValueSet, optionally providing the index (1-based) for each iteration in CountVar.";\

ProductForEach::usage="ProductForEach[Expr,{Var,ValueSet},...] or ProductForEach[Expr,{Var,CountVar,ValueSet},...]  returns a product of Expr, evaluated for each value of Var in the list ValueSet, optionally providing the index (1-based) for each iteration in CountVar.";



DoForEach::notlist="Value set must be a list (`1`).";
TableForEach::notlist="Value set must be a list (`1`).";



Begin["`Private`"];







SetAttributes[DoForEach,HoldAll];
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
            {Var},
            Var=ValueSet[[i]];
            Expr
            ],
          {i,1,Length[ValueSet]}
          ]
        ]
      
      ];
DoForEach[Expr_,{VarMultiplet:{_Symbol..},CountVar_Symbol,ValueSet_}]:=
    Module[
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
            VarMultiplet,
            VarMultiplet=ValueSet[[i]];
            Expr
            ],
          {i,1,Length[ValueSet]}
          ]
        ]
      
      ];
DoForEach[Expr_,{Var:({_Symbol..}|_Symbol),ValueSet_}]:=
  DoForEach[Expr,{Var,Dummy,ValueSet}];
DoForEach[Expr_,FirstIterator_List,RestIteratorSeq__List]:=
  DoForEach[DoForEach[Expr,RestIteratorSeq],FirstIterator];



Options[TableForEach]={Flatten\[Rule]False};
SetAttributes[TableForEach,HoldAll];
TableForEach[Expr_,{Var_Symbol,CountVar_Symbol,ValueSet_},Opts___?OptionQ]:=
    Module[
      {i},
      
      If[
        Head[ValueSet]=!=List,
        Message[TableForEach::notlist,ValueSet]
        ];
      
      Block[
        {CountVar},
        Table[
          CountVar=i;
          Block[
            {Var},
            Var=ValueSet[[i]];
            Expr
            ],
          {i,1,Length[ValueSet]}
          ]
        ]
      ];
TableForEach[Expr_,{VarMultiplet:{_Symbol..},CountVar_Symbol,ValueSet_},
      Opts___?OptionQ]:=Module[
      {i},
      
      If[
        Head[ValueSet]=!=List,
        Message[TableForEach::notlist,ValueSet]
        ];
      
      Block[
        {CountVar},
        Table[
          CountVar=i;
          Block[
            VarMultiplet,
            VarMultiplet=ValueSet[[i]];
            Expr
            ],
          {i,1,Length[ValueSet]}
          ]
        ]
      ];
TableForEach[Expr_,{Var:({_Symbol..}|_Symbol),ValueSet_},Opts___?OptionQ]:=
    TableForEach[Expr,{Var,Dummy,ValueSet},Opts];
TableForEach[
      Expr_,
      FirstIterator:{({_Symbol..}|_Symbol),__},
      RestIteratorSeq:({({_Symbol..}|_Symbol),__}..),
      Opts___?OptionQ
      ]:=Module[
      {FullOpts=Flatten[{Opts,Options[TableForEach]}]},
      If[Flatten/.FullOpts,
        Flatten[
          TableForEach[TableForEach[Expr,RestIteratorSeq,Opts],FirstIterator],
          1],
        TableForEach[TableForEach[Expr,RestIteratorSeq],FirstIterator]
        ]
      ];



SetAttributes[SumForEach,HoldAll];
SumForEach[Expr_,Iterators___]:=
    Plus@@TableForEach[Expr,Iterators,Flatten\[Rule]True];
SetAttributes[ProductForEach,HoldAll];
ProductForEach[Expr_,Iterators___]:=
    Times@@TableForEach[Expr,Iterators,Flatten\[Rule]True];



End[];



Protect[Evaluate[$Context<>"*"]];
EndPackage[];

