/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     ID = 258,
     INT_VAL = 259,
     FLOAT_VAL = 260,
     INFINITY = 261,
     COEFF = 262,
     COVER = 263,
     OBJ = 264,
     DEFAULT = 265,
     FROM = 266,
     TO = 267,
     TO_COME = 268,
     MODELTYPE = 269,
     NET_IN = 270,
     NET_OUT = 271,
     DIMEN = 272,
     ORDERED = 273,
     CIRCULAR = 274,
     REVERSED = 275,
     SYMBOLIC = 276,
     ARC = 277,
     INTEGER = 278,
     BINARY = 279,
     CHECK = 280,
     CLOSE = 281,
     DISPLAY = 282,
     DROP = 283,
     INCLUDE = 284,
     PRINT = 285,
     PRINTF = 286,
     QUIT = 287,
     RESET = 288,
     RESTORE = 289,
     SOLVE = 290,
     UPDATE = 291,
     WRITE = 292,
     SHELL = 293,
     MODEL = 294,
     DATA = 295,
     OPTION = 296,
     LET = 297,
     SOLUTION = 298,
     FIX = 299,
     UNFIX = 300,
     END = 301,
     FUNCTION = 302,
     PIPE = 303,
     FORMAT = 304,
     SETOF = 305,
     BY = 306,
     LESS = 307,
     MOD = 308,
     DIV = 309,
     MIN = 310,
     MAX = 311,
     IF = 312,
     THEN = 313,
     ELSE = 314,
     AND = 315,
     OR = 316,
     EXISTS = 317,
     FORALL = 318,
     NOT = 319,
     WITHIN = 320,
     WHILE = 321,
     REPEAT = 322,
     FOR = 323,
     CARD = 324,
     NEXT = 325,
     NEXTW = 326,
     PREV = 327,
     PREVW = 328,
     FIRST = 329,
     LAST = 330,
     MEMBER = 331,
     ORD = 332,
     ORD_ZERO = 333,
     VAR = 334,
     PARAM = 335,
     SET = 336,
     MAXIMIZE = 337,
     MINIMIZE = 338,
     OBJECTIVE = 339,
     SUBJECTTO = 340,
     SUM = 341,
     PROD = 342,
     IN = 343,
     POWER = 344,
     NE = 345,
     LE = 346,
     GE = 347,
     EQ = 348,
     LT = 349,
     GT = 350,
     UNION = 351,
     DIFF = 352,
     CROSS = 353,
     INTER = 354,
     SYMDIFF = 355,
     LBRACE = 356,
     RBRACE = 357,
     COMMA = 358,
     SEMICOLON = 359,
     LSBRACKET = 360,
     RSBRACKET = 361,
     COLON = 362,
     LBRACKET = 363,
     RBRACKET = 364,
     DEFINED = 365,
     LOGICAL_OR = 366,
     LOGICAL_AND = 367,
     ELLIPSE = 368,
     PUBLIC = 369,
     CORE = 370,
     DOT = 371,
     BEG = 372,
     TIMESTAGE = 373,
     RANDOM = 374,
     SUFFIX = 375,
     BLOCK = 376,
     IDREF = 377,
     IDREFM = 378,
     SBLOCK = 379,
     USING = 380,
     DETERMINISTIC = 381,
     EXPECTATION = 382,
     STOCHASTIC = 383,
     STAGES = 384,
     STAGE = 385,
     NODE = 386,
     TR = 387,
     ASSIGN = 388,
     TOKPARAMSPECLIST = 389,
     TOKPARAMTEMPLATE = 390,
     TOKVALUETABLELIST = 391,
     TOKVALUETABLE = 392,
     CHARACTER_STRING = 393,
     TOKSETSPEC = 394
   };
#endif
/* Tokens.  */
#define ID 258
#define INT_VAL 259
#define FLOAT_VAL 260
#define INFINITY 261
#define COEFF 262
#define COVER 263
#define OBJ 264
#define DEFAULT 265
#define FROM 266
#define TO 267
#define TO_COME 268
#define MODELTYPE 269
#define NET_IN 270
#define NET_OUT 271
#define DIMEN 272
#define ORDERED 273
#define CIRCULAR 274
#define REVERSED 275
#define SYMBOLIC 276
#define ARC 277
#define INTEGER 278
#define BINARY 279
#define CHECK 280
#define CLOSE 281
#define DISPLAY 282
#define DROP 283
#define INCLUDE 284
#define PRINT 285
#define PRINTF 286
#define QUIT 287
#define RESET 288
#define RESTORE 289
#define SOLVE 290
#define UPDATE 291
#define WRITE 292
#define SHELL 293
#define MODEL 294
#define DATA 295
#define OPTION 296
#define LET 297
#define SOLUTION 298
#define FIX 299
#define UNFIX 300
#define END 301
#define FUNCTION 302
#define PIPE 303
#define FORMAT 304
#define SETOF 305
#define BY 306
#define LESS 307
#define MOD 308
#define DIV 309
#define MIN 310
#define MAX 311
#define IF 312
#define THEN 313
#define ELSE 314
#define AND 315
#define OR 316
#define EXISTS 317
#define FORALL 318
#define NOT 319
#define WITHIN 320
#define WHILE 321
#define REPEAT 322
#define FOR 323
#define CARD 324
#define NEXT 325
#define NEXTW 326
#define PREV 327
#define PREVW 328
#define FIRST 329
#define LAST 330
#define MEMBER 331
#define ORD 332
#define ORD_ZERO 333
#define VAR 334
#define PARAM 335
#define SET 336
#define MAXIMIZE 337
#define MINIMIZE 338
#define OBJECTIVE 339
#define SUBJECTTO 340
#define SUM 341
#define PROD 342
#define IN 343
#define POWER 344
#define NE 345
#define LE 346
#define GE 347
#define EQ 348
#define LT 349
#define GT 350
#define UNION 351
#define DIFF 352
#define CROSS 353
#define INTER 354
#define SYMDIFF 355
#define LBRACE 356
#define RBRACE 357
#define COMMA 358
#define SEMICOLON 359
#define LSBRACKET 360
#define RSBRACKET 361
#define COLON 362
#define LBRACKET 363
#define RBRACKET 364
#define DEFINED 365
#define LOGICAL_OR 366
#define LOGICAL_AND 367
#define ELLIPSE 368
#define PUBLIC 369
#define CORE 370
#define DOT 371
#define BEG 372
#define TIMESTAGE 373
#define RANDOM 374
#define SUFFIX 375
#define BLOCK 376
#define IDREF 377
#define IDREFM 378
#define SBLOCK 379
#define USING 380
#define DETERMINISTIC 381
#define EXPECTATION 382
#define STOCHASTIC 383
#define STAGES 384
#define STAGE 385
#define NODE 386
#define TR 387
#define ASSIGN 388
#define TOKPARAMSPECLIST 389
#define TOKPARAMTEMPLATE 390
#define TOKVALUETABLELIST 391
#define TOKVALUETABLE 392
#define CHARACTER_STRING 393
#define TOKSETSPEC 394




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 43 "data.tab.ypp"
{
   long *ival;
   double *fval;
   char *string;
   SyntaxNode *opPtr;
   SyntaxNodeIx *opPtrIx;
}
/* Line 1489 of yacc.c.  */
#line 335 "data.tab.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE datalval;

