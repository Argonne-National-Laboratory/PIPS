/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse dataparse
#define yylex   datalex
#define yyerror dataerror
#define yylval  datalval
#define yychar  datachar
#define yydebug datadebug
#define yynerrs datanerrs


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




/* Copy the first part of user declarations.  */
#line 20 "data.tab.ypp"

   #include "AmplModel.h"
   #include "CompDescrParam.h"
   #include "DataNodes.h"
   #include "GlobalVariables.h"
   #include "ModelComp.h"
   #include "nodes.h"
   #include <cassert>
   #include <cstdlib>
   #include <iostream>

   #define YYERROR_VERBOSE

   using namespace std;

   int datalex(void);
   static void dataerror(const char *s);
   extern FILE *datain;
   extern int datalineno;
   static AmplModel *root;


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 1
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

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
/* Line 187 of yacc.c.  */
#line 412 "data.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 425 "data.tab.cpp"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   158

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  143
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  27
/* YYNRULES -- Number of rules.  */
#define YYNRULES  54
/* YYNRULES -- Number of states.  */
#define YYNSTATES  85

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   394

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,   140,   141,     2,   142,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    17,    18,    20,
      21,    24,    25,    29,    31,    35,    37,    39,    41,    43,
      46,    48,    51,    54,    57,    61,    64,    66,    68,    71,
      74,    78,    81,    83,    86,    88,    90,    92,    94,    96,
      98,   105,   113,   123,   124,   127,   130,   132,   135,   138,
     141,   142,   146,   148,   151
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     144,     0,    -1,    -1,   144,   145,    -1,   163,    -1,   146,
      -1,    81,     3,   147,   148,   104,    -1,    -1,   110,    -1,
      -1,   148,   154,    -1,    -1,   108,   150,   109,    -1,   151,
      -1,   150,   103,   151,    -1,   162,    -1,   140,    -1,   107,
      -1,   162,    -1,   152,   162,    -1,   155,    -1,   153,   155,
      -1,   149,   152,    -1,   149,   153,    -1,   156,   110,   157,
      -1,   132,   158,    -1,   158,    -1,   159,    -1,   157,   159,
      -1,   107,   152,    -1,   158,   107,   152,    -1,   152,   160,
      -1,   161,    -1,   160,   161,    -1,   141,    -1,   142,    -1,
       3,    -1,   138,    -1,     4,    -1,     5,    -1,    80,     3,
     164,   147,   165,   104,    -1,    80,   164,   107,   152,   110,
     152,   104,    -1,    80,   164,   107,   162,   107,   152,   110,
     152,   104,    -1,    -1,    10,     4,    -1,    10,     5,    -1,
     166,    -1,   165,   166,    -1,   167,   152,    -1,   167,   168,
      -1,    -1,   105,   150,   106,    -1,   169,    -1,   168,   169,
      -1,   156,   110,   152,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,    85,    85,    86,    89,    90,    94,   113,   114,   117,
     118,   124,   125,   128,   129,   132,   133,   134,   137,   138,
     141,   142,   145,   146,   149,   152,   156,   159,   160,   163,
     164,   170,   176,   177,   180,   181,   184,   185,   186,   187,
     190,   206,   208,   212,   213,   216,   221,   222,   225,   232,
     238,   239,   242,   243,   246
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ID", "INT_VAL", "FLOAT_VAL", "INFINITY",
  "COEFF", "COVER", "OBJ", "DEFAULT", "FROM", "TO", "TO_COME", "MODELTYPE",
  "NET_IN", "NET_OUT", "DIMEN", "ORDERED", "CIRCULAR", "REVERSED",
  "SYMBOLIC", "ARC", "INTEGER", "BINARY", "CHECK", "CLOSE", "DISPLAY",
  "DROP", "INCLUDE", "PRINT", "PRINTF", "QUIT", "RESET", "RESTORE",
  "SOLVE", "UPDATE", "WRITE", "SHELL", "MODEL", "DATA", "OPTION", "LET",
  "SOLUTION", "FIX", "UNFIX", "END", "FUNCTION", "PIPE", "FORMAT", "SETOF",
  "BY", "LESS", "MOD", "DIV", "MIN", "MAX", "IF", "THEN", "ELSE", "AND",
  "OR", "EXISTS", "FORALL", "NOT", "WITHIN", "WHILE", "REPEAT", "FOR",
  "CARD", "NEXT", "NEXTW", "PREV", "PREVW", "FIRST", "LAST", "MEMBER",
  "ORD", "ORD_ZERO", "VAR", "PARAM", "SET", "MAXIMIZE", "MINIMIZE",
  "OBJECTIVE", "SUBJECTTO", "SUM", "PROD", "IN", "POWER", "NE", "LE", "GE",
  "EQ", "LT", "GT", "UNION", "DIFF", "CROSS", "INTER", "SYMDIFF", "LBRACE",
  "RBRACE", "COMMA", "SEMICOLON", "LSBRACKET", "RSBRACKET", "COLON",
  "LBRACKET", "RBRACKET", "DEFINED", "LOGICAL_OR", "LOGICAL_AND",
  "ELLIPSE", "PUBLIC", "CORE", "DOT", "BEG", "TIMESTAGE", "RANDOM",
  "SUFFIX", "BLOCK", "IDREF", "IDREFM", "SBLOCK", "USING", "DETERMINISTIC",
  "EXPECTATION", "STOCHASTIC", "STAGES", "STAGE", "NODE", "TR", "ASSIGN",
  "TOKPARAMSPECLIST", "TOKPARAMTEMPLATE", "TOKVALUETABLELIST",
  "TOKVALUETABLE", "CHARACTER_STRING", "TOKSETSPEC", "'*'", "'+'", "'-'",
  "$accept", "statements", "statement", "setdef", "defined_opt",
  "setspec_list", "set_template_opt", "templ_item_list", "templ_item",
  "object_list", "member_table_list", "setspec", "member_table",
  "theader_row", "settablerow_list", "theader_list", "settablerow",
  "plusminuslist", "plusminus", "object", "paramdef", "paramdefault_opt",
  "paramspec_list", "paramspec", "paramtemplate_opt", "valuetable_list",
  "valuetable", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
      42,    43,    45
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   143,   144,   144,   145,   145,   146,   147,   147,   148,
     148,   149,   149,   150,   150,   151,   151,   151,   152,   152,
     153,   153,   154,   154,   155,   156,   156,   157,   157,   158,
     158,   159,   160,   160,   161,   161,   162,   162,   162,   162,
     163,   163,   163,   164,   164,   164,   165,   165,   166,   166,
     167,   167,   168,   168,   169
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     5,     0,     1,     0,
       2,     0,     3,     1,     3,     1,     1,     1,     1,     2,
       1,     2,     2,     2,     3,     2,     1,     1,     2,     2,
       3,     2,     1,     2,     1,     1,     1,     1,     1,     1,
       6,     7,     9,     0,     2,     2,     1,     2,     2,     2,
       0,     3,     1,     2,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    43,     0,     3,     5,     4,    43,     0,
       0,     7,     7,    44,    45,     0,     8,     9,    50,    36,
      38,    39,    37,     0,    18,    11,     0,    50,    46,     0,
       0,    19,     0,     6,     0,     0,    10,    17,    16,     0,
      13,    15,    40,    47,     0,     0,    48,     0,    26,    18,
      49,    52,     0,     0,     0,    22,    23,    20,     0,     0,
      51,    29,    25,     0,     0,    53,    41,     0,    12,    21,
       0,    14,    54,    30,     0,     0,    24,    27,    42,    34,
      35,    31,    32,    28,    33
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     5,     6,    17,    25,    35,    39,    40,    75,
      56,    36,    57,    47,    76,    48,    77,    81,    82,    31,
       7,    10,    27,    28,    29,    50,    51
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -104
static const yytype_int8 yypact[] =
{
    -104,     3,  -104,    23,    29,  -104,  -104,  -104,    30,    54,
     -64,   -65,   -65,  -104,  -104,    20,  -104,  -104,   -57,  -104,
    -104,  -104,  -104,     8,   -54,   -69,     2,   -44,  -104,     5,
      20,  -104,    20,  -104,     2,     5,  -104,  -104,  -104,   -72,
    -104,  -104,  -104,  -104,    20,   -52,    20,   -60,   -43,  -104,
    -103,  -104,    11,    14,   -73,    20,  -103,  -104,   -45,     2,
    -104,    20,   -43,    20,    20,  -104,  -104,    20,  -104,  -104,
      20,  -104,    20,    20,    17,    -3,    20,  -104,  -104,  -104,
    -104,   -79,  -104,  -104,  -104
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
    -104,  -104,  -104,  -104,    55,  -104,  -104,    34,    10,    22,
    -104,  -104,    16,    -7,  -104,    25,    -2,  -104,    -8,    12,
    -104,    69,  -104,    51,  -104,  -104,    31
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      19,    20,    21,     2,    44,    19,    20,    21,    19,    20,
      21,    19,    20,    21,    19,    20,    21,    19,    20,    21,
      19,    20,    21,    19,    20,    21,     8,    24,    58,    45,
      59,    59,    11,     9,    60,    33,    68,    23,    41,    34,
       9,    49,    49,    15,    49,    16,    41,    49,    26,    58,
      63,    46,    52,    32,    53,    44,    49,    55,    13,    14,
      42,    26,    79,    80,    64,    70,    61,    18,    54,    71,
      62,    41,    69,    84,    83,    49,    49,    12,    43,    49,
       0,    65,    49,     3,     4,    72,    73,     0,    49,    74,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    37,
       0,     0,    44,     0,     0,    66,     0,     0,    30,     0,
       0,    78,     0,     0,    67,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    22,     0,    45,    79,    80,
      22,     0,    38,    22,     0,     0,    22,     0,     0,    22,
       0,     0,    22,     0,     0,    22,     0,     0,    22
};

static const yytype_int16 yycheck[] =
{
       3,     4,     5,     0,   107,     3,     4,     5,     3,     4,
       5,     3,     4,     5,     3,     4,     5,     3,     4,     5,
       3,     4,     5,     3,     4,     5,     3,    15,    35,   132,
     103,   103,     3,    10,   106,   104,   109,    15,    26,   108,
      10,    29,    30,   107,    32,   110,    34,    35,   105,    56,
     110,    29,    30,   107,    32,   107,    44,    35,     4,     5,
     104,   105,   141,   142,   107,   110,    44,    12,    34,    59,
      45,    59,    56,    81,    76,    63,    64,     8,    27,    67,
      -1,    50,    70,    80,    81,    63,    64,    -1,    76,    67,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   107,
      -1,    -1,   107,    -1,    -1,   104,    -1,    -1,   110,    -1,
      -1,   104,    -1,    -1,   110,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   138,    -1,   132,   141,   142,
     138,    -1,   140,   138,    -1,    -1,   138,    -1,    -1,   138,
      -1,    -1,   138,    -1,    -1,   138,    -1,    -1,   138
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,   144,     0,    80,    81,   145,   146,   163,     3,    10,
     164,     3,   164,     4,     5,   107,   110,   147,   147,     3,
       4,     5,   138,   152,   162,   148,   105,   165,   166,   167,
     110,   162,   107,   104,   108,   149,   154,   107,   140,   150,
     151,   162,   104,   166,   107,   132,   152,   156,   158,   162,
     168,   169,   152,   152,   150,   152,   153,   155,   156,   103,
     106,   152,   158,   110,   107,   169,   104,   110,   109,   155,
     110,   151,   152,   152,   152,   152,   157,   159,   104,   141,
     142,   160,   161,   159,   161
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 6:
#line 94 "data.tab.ypp"
    {
            // Find symbol
            const SymbolTable::Entry *entry = root->symbol_table.findSymbol((yyvsp[(2) - (5)].string));
            if (!entry) {
               cerr << "Data given for set " << (yyvsp[(2) - (5)].string) << " which has not been"
                  "declared at top level!" << endl;
               exit(1);
            }
            if (!entry->isType(SymbolTable::ST_SET)) {
               cerr << "Set data given for variable " << (yyvsp[(2) - (5)].string) << " which was "
                  "not declared a set!" << endl;
               exit(1);
            }
            assert((yyvsp[(4) - (5)].opPtr)->nchild() == 1);
            SetSpec *tmp = (SetSpec*) (yyvsp[(4) - (5)].opPtr)->front();
            entry->mc->setValue(new Set(*(tmp->list)));
         }
    break;

  case 9:
#line 117 "data.tab.ypp"
    { (yyval.opPtr)=NULL; }
    break;

  case 10:
#line 118 "data.tab.ypp"
    { 
                  (yyval.opPtr) = addItemToListOrCreate(' ', (ListNode*)(yyvsp[(1) - (2)].opPtr), (yyvsp[(2) - (2)].opPtr)); 
               }
    break;

  case 11:
#line 124 "data.tab.ypp"
    { (yyval.opPtr)=NULL; }
    break;

  case 12:
#line 125 "data.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(2) - (3)].opPtr); }
    break;

  case 16:
#line 133 "data.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode('*'); }
    break;

  case 17:
#line 134 "data.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(COLON); }
    break;

  case 18:
#line 137 "data.tab.ypp"
    { (yyval.opPtr) = new ListNode(' ', (yyvsp[(1) - (1)].opPtr)); }
    break;

  case 19:
#line 138 "data.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(1) - (2)].opPtr)->push_back((yyvsp[(2) - (2)].opPtr)); }
    break;

  case 22:
#line 145 "data.tab.ypp"
    { (yyval.opPtr) = new SetSpec((yyvsp[(1) - (2)].opPtr), (ListNode*)(yyvsp[(2) - (2)].opPtr)); }
    break;

  case 25:
#line 152 "data.tab.ypp"
    {
                  cerr << "TR not yet supported." << endl;
                  (yyval.opPtr) = (yyvsp[(2) - (2)].opPtr);
               }
    break;

  case 29:
#line 163 "data.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(2) - (2)].opPtr); }
    break;

  case 30:
#line 164 "data.tab.ypp"
    {
                  cerr << "Multiple theaders not yet supported." << endl;
                  (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr);
               }
    break;

  case 31:
#line 170 "data.tab.ypp"
    {
                  cerr << "set table definitions not yet supported." << endl;
                  (yyval.opPtr) = NULL;
               }
    break;

  case 32:
#line 176 "data.tab.ypp"
    { (yyval.opPtr) = NULL; }
    break;

  case 33:
#line 177 "data.tab.ypp"
    { (yyval.opPtr) = NULL; }
    break;

  case 36:
#line 184 "data.tab.ypp"
    { (yyval.opPtr) = new IDNode((yyvsp[(1) - (1)].string)); }
    break;

  case 37:
#line 185 "data.tab.ypp"
    { (yyval.opPtr) = new IDNode((yyvsp[(1) - (1)].string)); }
    break;

  case 38:
#line 186 "data.tab.ypp"
    { (yyval.opPtr) = new ValueNode<long>(*(yyvsp[(1) - (1)].ival)); }
    break;

  case 39:
#line 187 "data.tab.ypp"
    { (yyval.opPtr) = new ValueNode<double>(*(yyvsp[(1) - (1)].fval)); }
    break;

  case 40:
#line 190 "data.tab.ypp"
    {
               const SymbolTable::Entry *entry =
                  root->symbol_table.findSymbol((yyvsp[(2) - (6)].string));
                if (entry==NULL){
                  cerr << "Data given for parameter " << (yyvsp[(2) - (6)].string) << " which has "
                    "not beendeclared at top level!" << endl;
                  exit(1);
                }
                if (!entry->isType(SymbolTable::ST_PARAM)) {
                  cerr << "Param data given for variable " << (yyvsp[(2) - (6)].string) << " which "
                    "was not delcared a param!" << endl;
                    exit(1);
                }
                // reference has been found
                entry->mc->setValue(new CompDescrParam(entry->mc, (yyvsp[(5) - (6)].opPtr)));
            }
    break;

  case 44:
#line 213 "data.tab.ypp"
    {
                     cerr << "Default values not yet supported." << endl;
                  }
    break;

  case 45:
#line 216 "data.tab.ypp"
    {
                     cerr << "Default values not yet supported." << endl;
                  }
    break;

  case 46:
#line 221 "data.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(TOKPARAMSPECLIST, (yyvsp[(1) - (1)].opPtr)); }
    break;

  case 47:
#line 222 "data.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(1) - (2)].opPtr)->push_back((yyvsp[(2) - (2)].opPtr)); }
    break;

  case 48:
#line 225 "data.tab.ypp"
    {
               if((yyvsp[(1) - (2)].opPtr)) {
                  (yyval.opPtr) = new SyntaxNode(TOKPARAMTEMPLATE, (yyvsp[(1) - (2)].opPtr), (yyvsp[(2) - (2)].opPtr));
               } else {
                  (yyval.opPtr) = (yyvsp[(2) - (2)].opPtr);
               }
            }
    break;

  case 50:
#line 238 "data.tab.ypp"
    { (yyval.opPtr)=NULL; }
    break;

  case 51:
#line 239 "data.tab.ypp"
    { (yyval.opPtr)=(yyvsp[(2) - (3)].opPtr); }
    break;

  case 52:
#line 242 "data.tab.ypp"
    { new SyntaxNode(TOKVALUETABLELIST, (yyvsp[(1) - (1)].opPtr)); }
    break;

  case 53:
#line 243 "data.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(1) - (2)].opPtr)->push_back((yyvsp[(2) - (2)].opPtr)); }
    break;

  case 54:
#line 246 "data.tab.ypp"
    {
               if((yyvsp[(3) - (3)].opPtr)->nchild() % ((yyvsp[(1) - (3)].opPtr)->nchild()+1) != 0) {
                  cerr << "DATA: error in line " << datalineno << ":" << endl;
                  cerr << "Length of value table (" << (yyvsp[(3) - (3)].opPtr)->nchild() << ")" <<
                     " not divisable by column labels +1 (" <<
                     (yyvsp[(1) - (3)].opPtr)->nchild()+1 << ")" << endl;
                  exit(1);
               }
               (yyval.opPtr) = new SyntaxNode(TOKVALUETABLE, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
            }
    break;


/* Line 1267 of yacc.c.  */
#line 1966 "data.tab.cpp"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 259 "data.tab.ypp"


void dataerror(const char *s) {
  cerr << "DATA: " << s << " on line " << datalineno << "\n";
  exit(1);
}

void parse_data(AmplModel *current_model, const string& datafilename) {
        
  if (GlobalVariables::prtLvl >= PRINT_LOG) {
    cout << "===============================================================\n";
    cout << " Start parsing data file: " << datafilename << "\n";
    cout << "===============================================================\n";
  }

  datain = fopen(datafilename.c_str(), "r");
  root = current_model;
  if (datain==NULL){
    cout << "ERROR: Data file '" << datafilename << "' not found\n";
    exit(1);
  }
          
  dataparse();
  if (GlobalVariables::prtLvl >= PRINT_LOG) {
    cout << "===============================================================\n";
    cout << " Finished parsing data file\n";
    cout << "===============================================================\n";
  }
}

