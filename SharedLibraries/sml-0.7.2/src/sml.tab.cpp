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
     DOTDOT = 270,
     NET_IN = 271,
     NET_OUT = 272,
     DIMEN = 273,
     ORDERED = 274,
     CIRCULAR = 275,
     REVERSED = 276,
     SYMBOLIC = 277,
     ARC = 278,
     INTEGER = 279,
     BINARY = 280,
     CHECK = 281,
     CLOSE = 282,
     DISPLAY = 283,
     DROP = 284,
     INCLUDE = 285,
     PRINT = 286,
     PRINTF = 287,
     QUIT = 288,
     RESET = 289,
     RESTORE = 290,
     SOLVE = 291,
     UPDATE = 292,
     WRITE = 293,
     SHELL = 294,
     MODEL = 295,
     DATA = 296,
     OPTION = 297,
     LET = 298,
     SOLUTION = 299,
     FIX = 300,
     UNFIX = 301,
     END = 302,
     FUNCTION = 303,
     PIPE = 304,
     FORMAT = 305,
     SETOF = 306,
     BY = 307,
     LESS = 308,
     MOD = 309,
     DIV = 310,
     MIN = 311,
     MAX = 312,
     IF = 313,
     THEN = 314,
     ELSE = 315,
     AND = 316,
     OR = 317,
     EXISTS = 318,
     FORALL = 319,
     NOT = 320,
     WITHIN = 321,
     WHILE = 322,
     REPEAT = 323,
     FOR = 324,
     CARD = 325,
     NEXT = 326,
     NEXTW = 327,
     PREV = 328,
     PREVW = 329,
     FIRST = 330,
     LAST = 331,
     MEMBER = 332,
     ORD = 333,
     ORD_ZERO = 334,
     VAR = 335,
     PARAM = 336,
     SET = 337,
     MAXIMIZE = 338,
     MINIMIZE = 339,
     OBJECTIVE = 340,
     SUBJECTTO = 341,
     SUM = 342,
     PROD = 343,
     IN = 344,
     POWER = 345,
     NE = 346,
     LE = 347,
     GE = 348,
     EQ = 349,
     LT = 350,
     GT = 351,
     UNION = 352,
     DIFF = 353,
     CROSS = 354,
     INTER = 355,
     SYMDIFF = 356,
     LBRACE = 357,
     RBRACE = 358,
     COMMA = 359,
     SEMICOLON = 360,
     LSBRACKET = 361,
     RSBRACKET = 362,
     COLON = 363,
     LBRACKET = 364,
     RBRACKET = 365,
     DEFINED = 366,
     LOGICAL_OR = 367,
     LOGICAL_AND = 368,
     ELLIPSE = 369,
     DOT = 370,
     SUFFIX = 371,
     BLOCK = 372,
     IDREF = 373,
     IDREFM = 374,
     STAGE = 375,
     NODE = 376,
     USING = 377,
     DETERMINISTIC = 378,
     EXPECTATION = 379,
     STOCHASTIC = 380,
     STAGES = 381,
     ANCESTOR = 382,
     ASSIGN = 383
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
#define DOTDOT 270
#define NET_IN 271
#define NET_OUT 272
#define DIMEN 273
#define ORDERED 274
#define CIRCULAR 275
#define REVERSED 276
#define SYMBOLIC 277
#define ARC 278
#define INTEGER 279
#define BINARY 280
#define CHECK 281
#define CLOSE 282
#define DISPLAY 283
#define DROP 284
#define INCLUDE 285
#define PRINT 286
#define PRINTF 287
#define QUIT 288
#define RESET 289
#define RESTORE 290
#define SOLVE 291
#define UPDATE 292
#define WRITE 293
#define SHELL 294
#define MODEL 295
#define DATA 296
#define OPTION 297
#define LET 298
#define SOLUTION 299
#define FIX 300
#define UNFIX 301
#define END 302
#define FUNCTION 303
#define PIPE 304
#define FORMAT 305
#define SETOF 306
#define BY 307
#define LESS 308
#define MOD 309
#define DIV 310
#define MIN 311
#define MAX 312
#define IF 313
#define THEN 314
#define ELSE 315
#define AND 316
#define OR 317
#define EXISTS 318
#define FORALL 319
#define NOT 320
#define WITHIN 321
#define WHILE 322
#define REPEAT 323
#define FOR 324
#define CARD 325
#define NEXT 326
#define NEXTW 327
#define PREV 328
#define PREVW 329
#define FIRST 330
#define LAST 331
#define MEMBER 332
#define ORD 333
#define ORD_ZERO 334
#define VAR 335
#define PARAM 336
#define SET 337
#define MAXIMIZE 338
#define MINIMIZE 339
#define OBJECTIVE 340
#define SUBJECTTO 341
#define SUM 342
#define PROD 343
#define IN 344
#define POWER 345
#define NE 346
#define LE 347
#define GE 348
#define EQ 349
#define LT 350
#define GT 351
#define UNION 352
#define DIFF 353
#define CROSS 354
#define INTER 355
#define SYMDIFF 356
#define LBRACE 357
#define RBRACE 358
#define COMMA 359
#define SEMICOLON 360
#define LSBRACKET 361
#define RSBRACKET 362
#define COLON 363
#define LBRACKET 364
#define RBRACKET 365
#define DEFINED 366
#define LOGICAL_OR 367
#define LOGICAL_AND 368
#define ELLIPSE 369
#define DOT 370
#define SUFFIX 371
#define BLOCK 372
#define IDREF 373
#define IDREFM 374
#define STAGE 375
#define NODE 376
#define USING 377
#define DETERMINISTIC 378
#define EXPECTATION 379
#define STOCHASTIC 380
#define STAGES 381
#define ANCESTOR 382
#define ASSIGN 383




/* Copy the first part of user declarations.  */
#line 19 "sml.tab.ypp"

   #define YYERROR_VERBOSE
   #include <stdio.h>
   #include <stdlib.h>
   #include <assert.h>
   #include <iostream>
   #include "nodes.h"
   #include "GlobalVariables.h"
   #include "StochModel.h"
   #include "StochModelComp.h"
   #include "symtab.h"
   #include "SetNode.h"
   #include "AmplModel.h"

   #ifdef HAVE_DIRECT_H
   #include <direct.h> // for chdir() under MinGW
   #endif

   using namespace std;

   void add_indexing(SyntaxNodeIx *indexing);
   void rem_indexing(SyntaxNodeIx *indexing);
   void begin_model(char *name, SyntaxNode *indexing);
   void end_model(char *name);
   void begin_smodel(char *name, SyntaxNode *indexing, SyntaxNode *stochsets);
   void end_smodel(char *name);

   extern int yylineno;
   int yylex(void);
   void yyerror(const char *s);

   static AmplModel *current_model;    /* this is the model currently active */
                                       /* this is the GLOBAL context */
   static AmplModel *local_context;    /* this is the LOCAL context */
   SyntaxNodeIx *list_of_indexing[20];    /* list of currently applicable 
                                                indexing expressions */
   int n_indexing;
   /* ---------------- stochastic global variables:------------------------ */
   static bool is_stoch_model;      /* true if inside stochastic model def */
   /* these are set by global stocastic modifier commands */
   static bool is_deterministic_glo;
   static SyntaxNode *stages_glo;
   extern FILE *yyin;

   void addStochInfo(ModelComp *newmc, SyntaxNode*);


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
#line 67 "sml.tab.ypp"
{
  int optype;
  long *ival;
  double *fval;
  char *string;
  SyntaxNode *opPtr;
  SyntaxNodeIx *opPtrIx;
}
/* Line 187 of yacc.c.  */
#line 408 "sml.tab.cpp"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 421 "sml.tab.cpp"

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
#define YYLAST   741

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  134
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  56
/* YYNRULES -- Number of rules.  */
#define YYNRULES  163
/* YYNRULES -- Number of states.  */
#define YYNSTATES  326

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   383

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   133,     2,     2,     2,     2,     2,     2,
       2,     2,   131,   129,     2,   130,     2,   132,     2,     2,
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
     125,   126,   127,   128
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    13,    17,    22,
      27,    31,    36,    40,    45,    55,    59,    64,    67,    72,
      76,    79,    82,    85,    88,    91,    93,    97,    99,   105,
     107,   109,   111,   113,   115,   117,   119,   123,   127,   133,
     137,   143,   147,   153,   157,   163,   167,   168,   176,   177,
     185,   186,   194,   195,   197,   200,   205,   208,   210,   214,
     216,   220,   224,   227,   230,   234,   236,   240,   242,   243,
     250,   251,   259,   260,   267,   268,   269,   271,   273,   276,
     279,   282,   285,   288,   290,   291,   293,   295,   299,   301,
     303,   305,   308,   311,   314,   317,   320,   321,   323,   326,
     327,   329,   331,   335,   337,   339,   341,   344,   347,   350,
     353,   356,   359,   363,   367,   369,   373,   377,   381,   385,
     389,   393,   395,   397,   399,   401,   403,   405,   409,   411,
     416,   423,   430,   432,   436,   438,   442,   446,   450,   453,
     457,   461,   465,   469,   473,   477,   480,   481,   486,   491,
     498,   503,   508,   513,   518,   520,   522,   524,   526,   528,
     530,   532,   534,   536
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     135,     0,    -1,    -1,   135,   146,    -1,   140,    -1,   137,
      -1,   144,    -1,   138,   135,   139,    -1,   138,   102,   135,
     103,    -1,   117,     3,   154,   108,    -1,    47,   117,   105,
      -1,    47,   117,     3,   105,    -1,   141,   135,   142,    -1,
     141,   102,   135,   103,    -1,   117,     3,   154,   125,   122,
     109,   157,   110,   108,    -1,    47,   143,   105,    -1,    47,
     143,     3,   105,    -1,   125,   117,    -1,   145,   102,   135,
     103,    -1,   126,   179,   108,    -1,   162,   105,    -1,   164,
     105,    -1,   166,   105,    -1,   150,   105,    -1,   147,   105,
      -1,   136,    -1,    30,     3,   105,    -1,    36,    -1,    43,
     156,   182,   111,   185,    -1,    92,    -1,    93,    -1,    95,
      -1,    96,    -1,   128,    -1,    94,    -1,    91,    -1,   185,
     128,   185,    -1,   185,    95,   185,    -1,   185,    95,   185,
      95,   185,    -1,   185,    92,   185,    -1,   185,    92,   185,
      92,   185,    -1,   185,    96,   185,    -1,   185,    96,   185,
      96,   185,    -1,   185,    93,   185,    -1,   185,    93,   185,
      93,   185,    -1,   185,    94,   185,    -1,    -1,    83,     3,
     154,   151,   175,   108,   185,    -1,    -1,    84,     3,   154,
     152,   175,   108,   185,    -1,    -1,    86,     3,   154,   153,
     175,   108,   149,    -1,    -1,   156,    -1,   102,   157,    -1,
     155,   108,   161,   103,    -1,   155,   103,    -1,   158,    -1,
     157,   104,   158,    -1,   179,    -1,     3,    89,   179,    -1,
     159,    89,   179,    -1,   160,   110,    -1,   109,     3,    -1,
     160,   104,     3,    -1,   185,    -1,   109,   161,   110,    -1,
     149,    -1,    -1,    82,     3,   154,   163,   175,   169,    -1,
      -1,   168,    81,     3,   154,   165,   175,   172,    -1,    -1,
      80,     3,   154,   167,   175,   176,    -1,    -1,    -1,   170,
      -1,   171,    -1,   170,   171,    -1,    18,     4,    -1,    66,
     179,    -1,   111,   179,    -1,    10,   179,    -1,    19,    -1,
      -1,   173,    -1,   174,    -1,   173,   104,   174,    -1,    25,
      -1,    24,    -1,    22,    -1,   148,   185,    -1,    89,   179,
      -1,   128,   185,    -1,    10,   185,    -1,   111,   185,    -1,
      -1,   123,    -1,   126,   179,    -1,    -1,   177,    -1,   178,
      -1,   177,   104,   178,    -1,    25,    -1,    24,    -1,    22,
      -1,    92,   185,    -1,    93,   185,    -1,   111,   185,    -1,
     128,   185,    -1,    10,   185,    -1,    89,   179,    -1,   116,
       3,   185,    -1,   102,   184,   103,    -1,   182,    -1,   179,
     180,   179,    -1,   179,   181,   179,    -1,   181,   156,   179,
      -1,   185,    15,   185,    -1,    51,   156,   182,    -1,   109,
     179,   110,    -1,    98,    -1,   101,    -1,    99,    -1,    97,
      -1,   100,    -1,   183,    -1,   182,   115,   183,    -1,     3,
      -1,     3,   106,   184,   107,    -1,     3,   109,     4,   105,
     184,   110,    -1,   127,   109,     4,   110,   115,   183,    -1,
     185,    -1,   184,   104,   185,    -1,   189,    -1,   109,   184,
     110,    -1,   185,   129,   185,    -1,   185,   130,   185,    -1,
     130,   185,    -1,   185,   131,   185,    -1,   185,   132,   185,
      -1,   185,    90,   185,    -1,   185,   114,   185,    -1,   185,
     112,   185,    -1,   185,   113,   185,    -1,   133,   185,    -1,
      -1,   188,   156,   186,   185,    -1,    58,   161,    59,   185,
      -1,    58,   161,    59,   185,    60,   185,    -1,    75,   109,
     179,   110,    -1,    76,   109,   179,   110,    -1,   124,   109,
     184,   110,    -1,   187,   109,   184,   110,    -1,    78,    -1,
      70,    -1,    87,    -1,    57,    -1,    56,    -1,    88,    -1,
       4,    -1,     5,    -1,   182,    -1,     6,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   120,   120,   121,   124,   125,   126,   154,   155,   161,
     167,   171,   206,   207,   216,   227,   231,   237,   245,   253,
     268,   269,   270,   271,   272,   273,   274,   282,   283,   286,
     287,   288,   289,   290,   291,   292,   295,   296,   297,   301,
     302,   306,   307,   311,   312,   316,   319,   319,   336,   336,
     353,   353,   371,   372,   375,   383,   389,   399,   400,   408,
     409,   412,   417,   420,   423,   428,   429,   431,   434,   434,
     457,   457,   476,   476,   497,   504,   505,   508,   509,   525,
     528,   531,   534,   537,   540,   541,   544,   548,   561,   562,
     563,   564,   565,   566,   567,   568,   571,   572,   581,   592,
     593,   597,   599,   605,   606,   607,   608,   609,   610,   611,
     612,   616,   617,   621,   622,   623,   626,   629,   633,   636,
     641,   649,   650,   651,   654,   655,   664,   668,   721,   724,
     730,   739,   764,   767,   780,   781,   782,   783,   784,   785,
     786,   787,   788,   789,   790,   791,   792,   792,   798,   799,
     800,   801,   802,   806,   811,   812,   815,   816,   817,   818,
     821,   824,   827,   830
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "ID", "INT_VAL", "FLOAT_VAL", "INFINITY",
  "COEFF", "COVER", "OBJ", "DEFAULT", "FROM", "TO", "TO_COME", "MODELTYPE",
  "DOTDOT", "NET_IN", "NET_OUT", "DIMEN", "ORDERED", "CIRCULAR",
  "REVERSED", "SYMBOLIC", "ARC", "INTEGER", "BINARY", "CHECK", "CLOSE",
  "DISPLAY", "DROP", "INCLUDE", "PRINT", "PRINTF", "QUIT", "RESET",
  "RESTORE", "SOLVE", "UPDATE", "WRITE", "SHELL", "MODEL", "DATA",
  "OPTION", "LET", "SOLUTION", "FIX", "UNFIX", "END", "FUNCTION", "PIPE",
  "FORMAT", "SETOF", "BY", "LESS", "MOD", "DIV", "MIN", "MAX", "IF",
  "THEN", "ELSE", "AND", "OR", "EXISTS", "FORALL", "NOT", "WITHIN",
  "WHILE", "REPEAT", "FOR", "CARD", "NEXT", "NEXTW", "PREV", "PREVW",
  "FIRST", "LAST", "MEMBER", "ORD", "ORD_ZERO", "VAR", "PARAM", "SET",
  "MAXIMIZE", "MINIMIZE", "OBJECTIVE", "SUBJECTTO", "SUM", "PROD", "IN",
  "POWER", "NE", "LE", "GE", "EQ", "LT", "GT", "UNION", "DIFF", "CROSS",
  "INTER", "SYMDIFF", "LBRACE", "RBRACE", "COMMA", "SEMICOLON",
  "LSBRACKET", "RSBRACKET", "COLON", "LBRACKET", "RBRACKET", "DEFINED",
  "LOGICAL_OR", "LOGICAL_AND", "ELLIPSE", "DOT", "SUFFIX", "BLOCK",
  "IDREF", "IDREFM", "STAGE", "NODE", "USING", "DETERMINISTIC",
  "EXPECTATION", "STOCHASTIC", "STAGES", "ANCESTOR", "ASSIGN", "'+'",
  "'-'", "'*'", "'/'", "'!'", "$accept", "statements", "block",
  "blockblock", "blockblockbegin", "blockblockend", "stochblock",
  "stochblockbegin", "stochblockend", "sblock_alias", "stageblock",
  "stageblock_start", "statement", "command", "relop", "equation", "cnstr",
  "@1", "@2", "@3", "indexing_opt", "start_indexing", "indexing",
  "setexpr_list", "setexpr_item", "mdim_dummy", "mdim_dummy_start",
  "lexpr", "setdef", "@4", "paramdef", "@5", "vardef", "@6", "qualifier",
  "setattributes_opt", "setattributes", "setattribute",
  "paramattributes_opt", "paramattributes", "paramattribute",
  "stochattr_opt", "varattributes_opt", "varattributes", "varattribute",
  "setexpression", "bsetop", "ubsetop", "identifier", "iditem",
  "expr_list", "expr", "@7", "func", "reduction_op", "value", 0
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
     375,   376,   377,   378,   379,   380,   381,   382,   383,    43,
      45,    42,    47,    33
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,   134,   135,   135,   136,   136,   136,   137,   137,   138,
     139,   139,   140,   140,   141,   142,   142,   143,   144,   145,
     146,   146,   146,   146,   146,   146,   146,   147,   147,   148,
     148,   148,   148,   148,   148,   148,   149,   149,   149,   149,
     149,   149,   149,   149,   149,   149,   151,   150,   152,   150,
     153,   150,   154,   154,   155,   156,   156,   157,   157,   158,
     158,   158,   159,   160,   160,   161,   161,   161,   163,   162,
     165,   164,   167,   166,   168,   169,   169,   170,   170,   171,
     171,   171,   171,   171,   172,   172,   173,   173,   174,   174,
     174,   174,   174,   174,   174,   174,   175,   175,   175,   176,
     176,   177,   177,   178,   178,   178,   178,   178,   178,   178,
     178,   178,   178,   179,   179,   179,   179,   179,   179,   179,
     179,   180,   180,   180,   181,   181,   182,   182,   183,   183,
     183,   183,   184,   184,   185,   185,   185,   185,   185,   185,
     185,   185,   185,   185,   185,   185,   186,   185,   185,   185,
     185,   185,   185,   185,   187,   187,   188,   188,   188,   188,
     189,   189,   189,   189
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     1,     3,     4,     4,
       3,     4,     3,     4,     9,     3,     4,     2,     4,     3,
       2,     2,     2,     2,     2,     1,     3,     1,     5,     1,
       1,     1,     1,     1,     1,     1,     3,     3,     5,     3,
       5,     3,     5,     3,     5,     3,     0,     7,     0,     7,
       0,     7,     0,     1,     2,     4,     2,     1,     3,     1,
       3,     3,     2,     2,     3,     1,     3,     1,     0,     6,
       0,     7,     0,     6,     0,     0,     1,     1,     2,     2,
       2,     2,     2,     1,     0,     1,     1,     3,     1,     1,
       1,     2,     2,     2,     2,     2,     0,     1,     2,     0,
       1,     1,     3,     1,     1,     1,     2,     2,     2,     2,
       2,     2,     3,     3,     1,     3,     3,     3,     3,     3,
       3,     1,     1,     1,     1,     1,     1,     3,     1,     4,
       6,     6,     1,     3,     1,     3,     3,     3,     2,     3,
       3,     3,     3,     3,     3,     2,     0,     4,     4,     6,
       4,     4,     4,     4,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,    74,     1,     0,    27,     0,     0,     0,     0,     0,
       0,     0,     0,    25,     5,     2,     4,     2,     6,     0,
       3,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    52,    52,    52,    52,    52,    52,   128,   160,   161,
     163,     0,   158,   157,     0,   155,     0,     0,   154,   156,
     159,   124,   125,     0,     0,     0,     0,     0,     0,     0,
       0,   114,   126,     0,     0,     0,   134,     2,    74,     2,
      74,     2,    24,    23,    20,    21,    22,     0,    26,   128,
       0,    54,    57,     0,     0,    59,    56,     0,     0,    72,
      53,    68,    46,    48,    50,     0,     0,     0,     0,     0,
      67,     0,   162,    65,     0,     0,     0,     0,   132,     0,
       0,   132,     0,     0,   138,   145,   121,   123,   122,    19,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   146,    74,     0,     7,    74,     0,
      12,    74,    52,     0,   128,     0,     0,     0,    62,     0,
       0,    96,    96,    96,    96,    96,     9,     0,     0,     0,
     119,     0,    65,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   113,     0,   120,   135,     0,     0,   115,   116,
     117,   127,   118,   141,   143,   144,   142,   136,   137,   139,
     140,     0,     0,     8,     0,    13,     0,     0,    18,    70,
      60,    58,    61,    64,    55,    28,    97,     0,    99,    75,
       0,     0,     0,     0,   129,     0,    66,   148,    39,    43,
      45,    37,    41,    36,   150,   151,   133,   152,     0,   153,
     147,     0,    10,    17,     0,    15,    96,    98,     0,   105,
     104,   103,     0,     0,     0,     0,     0,     0,    73,   100,
     101,     0,     0,    83,     0,     0,    69,    76,    77,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      11,    16,    84,   110,   111,   106,   107,   108,     0,   109,
       0,    82,    79,    80,    81,    78,    47,    49,    51,     0,
       0,   130,   149,    40,    44,    38,    42,   131,     0,    90,
      89,    88,     0,    35,    29,    30,    34,    31,    32,     0,
       0,     0,    71,    85,    86,   112,   102,     0,    94,    92,
      95,    93,    91,     0,    14,    87
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,    13,    14,    15,   137,    16,    17,   140,   197,
      18,    19,    20,    21,   311,   100,    22,   153,   154,   155,
      89,    29,    90,    81,    82,    83,    84,   101,    23,   152,
      24,   236,    25,   151,    26,   256,   257,   258,   312,   313,
     314,   208,   248,   249,   250,    85,   120,    60,   102,    62,
     110,    63,   192,    64,    65,    66
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -133
static const yytype_int16 yypact[] =
{
    -133,    40,  -133,     3,  -133,   -94,    47,    63,    71,    82,
      92,   102,   316,  -133,  -133,   -27,  -133,    -8,  -133,     4,
    -133,    24,    25,    32,    36,    41,    26,    42,   320,   -74,
       0,   -94,   -94,   -94,   -94,   -94,   -94,    50,  -133,  -133,
    -133,   -94,  -133,  -133,   429,  -133,    51,    58,  -133,  -133,
    -133,  -133,  -133,   467,   316,    59,    60,   467,   467,    54,
     -94,   -13,  -133,   173,    66,   -94,  -133,  -133,   466,  -133,
     484,  -133,  -133,  -133,  -133,  -133,  -133,   174,  -133,   -70,
     406,    76,  -133,    94,   -49,   199,  -133,   429,   -84,  -133,
    -133,  -133,  -133,  -133,  -133,   -80,   467,   177,     0,   429,
    -133,   127,    74,   251,   316,   316,   467,   -39,   609,   368,
     -41,   173,   467,   187,   -60,    78,  -133,  -133,  -133,  -133,
     316,   316,   316,     0,   467,   467,   467,   467,   467,   467,
     467,   467,   467,   467,  -133,   504,    81,  -133,   542,    68,
    -133,   566,   -94,   316,   -58,   320,   316,   196,  -133,    98,
     467,    35,    35,    35,    35,    35,  -133,    80,    72,    99,
      74,    95,   543,   467,   467,   467,   467,   467,   467,   467,
     391,   517,  -133,   467,  -133,  -133,   -26,    96,   199,   199,
     199,  -133,   609,   -65,   201,   -73,  -133,   -60,   -60,   -65,
     -65,   -22,   467,  -133,     6,  -133,   108,     7,  -133,  -133,
     199,  -133,   199,  -133,  -133,   609,  -133,   316,   160,    14,
     100,   119,   120,   117,  -133,   467,  -133,   180,   550,   577,
     609,   598,   588,   609,  -133,  -133,   609,  -133,   121,  -133,
     609,   128,  -133,  -133,   129,  -133,    35,   199,   467,  -133,
    -133,  -133,   316,   467,   467,   467,   229,   467,  -133,   133,
    -133,   316,   234,  -133,   316,   316,  -133,    14,  -133,   467,
     467,   467,   320,   -14,   467,   467,   467,   467,   467,     0,
    -133,  -133,   118,   609,   199,   609,   609,   609,   467,   609,
     160,   199,  -133,   199,   199,  -133,   609,   609,  -133,   251,
      28,  -133,   609,   609,   609,   609,   609,  -133,   467,  -133,
    -133,  -133,   316,  -133,  -133,  -133,  -133,  -133,  -133,   467,
     467,   467,  -133,   141,  -133,   609,  -133,   143,   609,   199,
     609,   609,   609,   118,  -133,  -133
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -133,    20,  -133,  -133,  -133,  -133,  -133,  -133,  -133,  -133,
    -133,  -133,  -133,  -133,  -133,     1,  -133,  -133,  -133,  -133,
     -21,  -133,     2,    -1,   122,  -133,  -133,   -61,  -133,  -133,
    -133,  -133,  -133,  -133,  -133,  -133,  -133,     8,  -133,  -133,
     -57,  -132,  -133,  -133,   -11,    -7,  -133,    -6,   -12,  -119,
     -52,    91,  -133,  -133,  -133,  -133
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -163
static const yytype_int16 yytable[] =
{
      61,   107,  -162,    37,   181,    59,    27,    30,    28,   231,
     234,    91,    92,    93,    94,    95,    61,   125,    88,   143,
     209,   210,   211,   212,   251,   125,   149,   150,   156,    86,
     125,   123,   252,   253,    87,    68,    96,    70,   161,    97,
       2,   128,    61,    98,   158,   157,   -63,   109,    96,   128,
      31,    97,   -63,   121,   128,   147,   129,   130,   131,   132,
     176,   148,   122,   173,   172,   173,    32,   134,    61,   175,
       3,   131,   132,   109,    33,    67,     4,  -162,   173,   121,
     254,   191,   173,     5,   227,    34,   160,   135,   229,   138,
     173,   141,    61,    61,    69,    35,   291,   170,   171,  -162,
    -162,  -162,   123,   121,   272,    36,    71,    77,    61,    61,
      61,   232,   235,   178,   179,   180,  -162,  -162,  -162,  -162,
       6,   199,     7,     8,     9,   255,    10,    56,   298,    72,
      73,    61,   145,    61,    61,   103,   200,    74,   317,   202,
     299,    75,   300,   301,   108,   111,    76,    78,   114,   115,
     297,    51,   116,   117,    52,   118,    96,    11,   206,    97,
     104,   207,   119,   263,   121,   121,    12,   105,   112,   113,
     238,   111,   121,   121,   121,   133,   173,   142,   103,   214,
     145,   159,   239,   146,   240,   241,   163,   108,   124,   123,
     162,   177,   128,   196,   121,    61,   121,   108,   194,   203,
     237,   204,   213,   108,   215,   216,   228,   302,   259,   303,
     304,   305,   306,   307,   308,   182,   183,   184,   185,   186,
     187,   188,   189,   190,   108,   233,   262,   260,   261,   309,
      61,   121,   278,   270,   271,   274,   269,   280,   282,    61,
     264,   205,    61,    61,   281,   323,   310,   283,   284,   242,
      61,   324,   243,   244,   217,   218,   219,   220,   221,   222,
     223,   290,   288,   125,   226,   285,   325,   201,   121,   316,
     125,   245,     0,     0,     0,   121,   246,   121,   121,     0,
       0,     0,     0,   230,     0,   126,   127,   128,   247,     0,
      61,   125,   126,   127,   128,   319,    51,   116,   117,    52,
     118,     0,   129,   130,   131,   132,   108,     0,     0,   129,
     130,   131,   132,   121,   127,   128,     0,     0,     0,    37,
      38,    39,    40,    79,    38,    39,    40,     0,     0,   273,
     129,   130,   131,   132,   275,   276,   277,     0,   279,     0,
       0,   125,     0,   164,   165,   166,   167,   168,     0,     0,
     286,   287,   289,     0,     0,   292,   293,   294,   295,   296,
       0,     0,     0,   126,   127,   128,     0,    41,     0,   315,
       0,    41,    42,    43,    44,     0,    42,    43,    44,   169,
     129,   130,   131,   132,     0,     0,    45,     0,     0,   318,
      45,    46,    47,     0,    48,    46,    47,     0,    48,     0,
     320,   321,   322,    49,    50,     0,     0,    49,    50,   144,
      38,    39,    40,    51,     0,     0,    52,    51,    53,     0,
      52,     0,    53,     0,     0,    54,     0,     0,     0,    80,
       0,     0,    37,    38,    39,    40,     0,     0,     0,     0,
      55,     0,     0,    56,    55,     0,    57,    56,     0,    58,
      57,     0,     0,    58,     0,     0,     0,    41,     0,     0,
       0,     0,    42,    43,    44,    51,   116,   117,    52,   118,
      37,    38,    39,    40,     0,     0,    45,     0,   174,     0,
       0,    46,    47,     0,    48,    42,    43,    44,    51,   116,
     117,    52,   118,    49,    50,     0,     3,     0,     0,    45,
       0,   224,     4,    51,    46,    47,    52,    48,    53,     5,
       0,     0,     0,   136,     3,    54,    49,    50,     0,     0,
       4,     0,     0,    42,    43,    44,     0,     5,     0,     0,
      55,   139,     0,    56,     3,     0,    57,    45,    99,    58,
       4,     0,    46,    47,     0,    48,     6,     5,     7,     8,
       9,     0,    10,    55,    49,    50,    56,     0,     0,    57,
       0,     0,    58,     0,     6,     0,     7,     8,     9,     0,
      10,     0,     3,     0,     0,     0,   106,     0,     4,     0,
       0,     0,     0,    11,     6,     5,     7,     8,     9,     0,
      10,    55,    12,     0,    56,     0,     3,    57,     0,     0,
      58,    11,     4,     0,     0,     0,     0,   193,     0,     5,
      12,     0,     0,     0,    51,   116,   117,    52,   118,     0,
       0,    11,     6,     0,     7,     8,     9,   225,    10,     0,
      12,     0,     0,   125,     0,   164,   165,   166,   167,   168,
     125,     0,   265,     0,     0,   195,     6,  -132,     7,     8,
       9,     0,    10,     0,     0,   126,   127,   128,     0,    11,
       0,     0,   126,   127,   128,     0,     0,   125,    12,   198,
     266,   169,   129,   130,   131,   132,     0,     0,   125,   129,
     130,   131,   132,    11,   268,     0,     0,     0,   125,   126,
     127,   128,    12,   267,     0,     0,     0,     0,     0,   125,
     126,   127,   128,     0,     0,     0,   129,   130,   131,   132,
     126,   127,   128,     0,     0,     0,     0,   129,   130,   131,
     132,   126,   127,   128,     0,     0,     0,   129,   130,   131,
     132,     0,     0,     0,     0,     0,     0,     0,   129,   130,
     131,   132
};

static const yytype_int16 yycheck[] =
{
      12,    53,    15,     3,   123,    12,     3,     5,   102,     3,
       3,    32,    33,    34,    35,    36,    28,    90,    30,    89,
     152,   153,   154,   155,    10,    90,    87,   111,   108,   103,
      90,   115,    18,    19,   108,    15,   106,    17,    99,   109,
       0,   114,    54,    41,    96,   125,   104,    54,   106,   114,
       3,   109,   110,    59,   114,   104,   129,   130,   131,   132,
     112,   110,    60,   104,   103,   104,     3,    65,    80,   110,
      30,   131,   132,    80,     3,   102,    36,    90,   104,    85,
      66,   133,   104,    43,   110,     3,    98,    67,   110,    69,
     104,    71,   104,   105,   102,     3,   110,   104,   105,   112,
     113,   114,   115,   109,   236,     3,   102,    81,   120,   121,
     122,   105,   105,   120,   121,   122,   129,   130,   131,   132,
      80,   142,    82,    83,    84,   111,    86,   127,    10,   105,
     105,   143,   104,   145,   146,    44,   143,   105,   110,   146,
      22,   105,    24,    25,    53,    54,   105,   105,    57,    58,
     269,    97,    98,    99,   100,   101,   106,   117,   123,   109,
     109,   126,   108,   215,   170,   171,   126,   109,   109,   109,
      10,    80,   178,   179,   180,   109,   104,     3,    87,   107,
     104,     4,    22,    89,    24,    25,    59,    96,    15,   115,
      99,     4,   114,   125,   200,   207,   202,   106,   117,     3,
     207,   103,   122,   112,   105,   110,   110,    89,   108,    91,
      92,    93,    94,    95,    96,   124,   125,   126,   127,   128,
     129,   130,   131,   132,   133,   117,   109,   108,   108,   111,
     242,   237,     3,   105,   105,   242,   115,   104,     4,   251,
      60,   150,   254,   255,   251,   104,   128,   254,   255,    89,
     262,   108,    92,    93,   163,   164,   165,   166,   167,   168,
     169,   262,   261,    90,   173,   257,   323,   145,   274,   280,
      90,   111,    -1,    -1,    -1,   281,   116,   283,   284,    -1,
      -1,    -1,    -1,   192,    -1,   112,   113,   114,   128,    -1,
     302,    90,   112,   113,   114,   302,    97,    98,    99,   100,
     101,    -1,   129,   130,   131,   132,   215,    -1,    -1,   129,
     130,   131,   132,   319,   113,   114,    -1,    -1,    -1,     3,
       4,     5,     6,     3,     4,     5,     6,    -1,    -1,   238,
     129,   130,   131,   132,   243,   244,   245,    -1,   247,    -1,
      -1,    90,    -1,    92,    93,    94,    95,    96,    -1,    -1,
     259,   260,   261,    -1,    -1,   264,   265,   266,   267,   268,
      -1,    -1,    -1,   112,   113,   114,    -1,    51,    -1,   278,
      -1,    51,    56,    57,    58,    -1,    56,    57,    58,   128,
     129,   130,   131,   132,    -1,    -1,    70,    -1,    -1,   298,
      70,    75,    76,    -1,    78,    75,    76,    -1,    78,    -1,
     309,   310,   311,    87,    88,    -1,    -1,    87,    88,     3,
       4,     5,     6,    97,    -1,    -1,   100,    97,   102,    -1,
     100,    -1,   102,    -1,    -1,   109,    -1,    -1,    -1,   109,
      -1,    -1,     3,     4,     5,     6,    -1,    -1,    -1,    -1,
     124,    -1,    -1,   127,   124,    -1,   130,   127,    -1,   133,
     130,    -1,    -1,   133,    -1,    -1,    -1,    51,    -1,    -1,
      -1,    -1,    56,    57,    58,    97,    98,    99,   100,   101,
       3,     4,     5,     6,    -1,    -1,    70,    -1,   110,    -1,
      -1,    75,    76,    -1,    78,    56,    57,    58,    97,    98,
      99,   100,   101,    87,    88,    -1,    30,    -1,    -1,    70,
      -1,   110,    36,    97,    75,    76,   100,    78,   102,    43,
      -1,    -1,    -1,    47,    30,   109,    87,    88,    -1,    -1,
      36,    -1,    -1,    56,    57,    58,    -1,    43,    -1,    -1,
     124,    47,    -1,   127,    30,    -1,   130,    70,   109,   133,
      36,    -1,    75,    76,    -1,    78,    80,    43,    82,    83,
      84,    -1,    86,   124,    87,    88,   127,    -1,    -1,   130,
      -1,    -1,   133,    -1,    80,    -1,    82,    83,    84,    -1,
      86,    -1,    30,    -1,    -1,    -1,   109,    -1,    36,    -1,
      -1,    -1,    -1,   117,    80,    43,    82,    83,    84,    -1,
      86,   124,   126,    -1,   127,    -1,    30,   130,    -1,    -1,
     133,   117,    36,    -1,    -1,    -1,    -1,   103,    -1,    43,
     126,    -1,    -1,    -1,    97,    98,    99,   100,   101,    -1,
      -1,   117,    80,    -1,    82,    83,    84,   110,    86,    -1,
     126,    -1,    -1,    90,    -1,    92,    93,    94,    95,    96,
      90,    -1,    92,    -1,    -1,   103,    80,   104,    82,    83,
      84,    -1,    86,    -1,    -1,   112,   113,   114,    -1,   117,
      -1,    -1,   112,   113,   114,    -1,    -1,    90,   126,   103,
      93,   128,   129,   130,   131,   132,    -1,    -1,    90,   129,
     130,   131,   132,   117,    96,    -1,    -1,    -1,    90,   112,
     113,   114,   126,    95,    -1,    -1,    -1,    -1,    -1,    90,
     112,   113,   114,    -1,    -1,    -1,   129,   130,   131,   132,
     112,   113,   114,    -1,    -1,    -1,    -1,   129,   130,   131,
     132,   112,   113,   114,    -1,    -1,    -1,   129,   130,   131,
     132,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   129,   130,
     131,   132
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,   135,     0,    30,    36,    43,    80,    82,    83,    84,
      86,   117,   126,   136,   137,   138,   140,   141,   144,   145,
     146,   147,   150,   162,   164,   166,   168,     3,   102,   155,
     156,     3,     3,     3,     3,     3,     3,     3,     4,     5,
       6,    51,    56,    57,    58,    70,    75,    76,    78,    87,
      88,    97,   100,   102,   109,   124,   127,   130,   133,   179,
     181,   182,   183,   185,   187,   188,   189,   102,   135,   102,
     135,   102,   105,   105,   105,   105,   105,    81,   105,     3,
     109,   157,   158,   159,   160,   179,   103,   108,   182,   154,
     156,   154,   154,   154,   154,   154,   106,   109,   156,   109,
     149,   161,   182,   185,   109,   109,   109,   184,   185,   179,
     184,   185,   109,   109,   185,   185,    98,    99,   101,   108,
     180,   181,   156,   115,    15,    90,   112,   113,   114,   129,
     130,   131,   132,   109,   156,   135,    47,   139,   135,    47,
     142,   135,     3,    89,     3,   104,    89,   104,   110,   161,
     111,   167,   163,   151,   152,   153,   108,   125,   184,     4,
     182,   161,   185,    59,    92,    93,    94,    95,    96,   128,
     179,   179,   103,   104,   110,   110,   184,     4,   179,   179,
     179,   183,   185,   185,   185,   185,   185,   185,   185,   185,
     185,   184,   186,   103,   117,   103,   125,   143,   103,   154,
     179,   158,   179,     3,   103,   185,   123,   126,   175,   175,
     175,   175,   175,   122,   107,   105,   110,   185,   185,   185,
     185,   185,   185,   185,   110,   110,   185,   110,   110,   110,
     185,     3,   105,   117,     3,   105,   165,   179,    10,    22,
      24,    25,    89,    92,    93,   111,   116,   128,   176,   177,
     178,    10,    18,    19,    66,   111,   169,   170,   171,   108,
     108,   108,   109,   184,    60,    92,    93,    95,    96,   115,
     105,   105,   175,   185,   179,   185,   185,   185,     3,   185,
     104,   179,     4,   179,   179,   171,   185,   185,   149,   185,
     157,   110,   185,   185,   185,   185,   185,   183,    10,    22,
      24,    25,    89,    91,    92,    93,    94,    95,    96,   111,
     128,   148,   172,   173,   174,   185,   178,   110,   185,   179,
     185,   185,   185,   104,   108,   174
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
        case 8:
#line 155 "sml.tab.ypp"
    {
               end_model(NULL);
               rem_indexing(NULL);
            }
    break;

  case 9:
#line 161 "sml.tab.ypp"
    {
                     add_indexing((yyvsp[(3) - (4)].opPtrIx));
                     begin_model((char*)(yyvsp[(2) - (4)].string), (yyvsp[(3) - (4)].opPtrIx));
                  }
    break;

  case 10:
#line 167 "sml.tab.ypp"
    {
                  end_model(NULL);
                  rem_indexing(NULL);
               }
    break;

  case 11:
#line 171 "sml.tab.ypp"
    {
                  end_model((char*)(yyvsp[(3) - (4)].string));
                  rem_indexing(NULL);
               }
    break;

  case 13:
#line 207 "sml.tab.ypp"
    {
               end_smodel(NULL);
               rem_indexing(NULL);
            }
    break;

  case 14:
#line 217 "sml.tab.ypp"
    {
                     add_indexing((yyvsp[(3) - (9)].opPtrIx));
                     begin_smodel((char*)(yyvsp[(2) - (9)].string), (yyvsp[(3) - (9)].opPtrIx), (yyvsp[(7) - (9)].opPtr));
                  }
    break;

  case 15:
#line 227 "sml.tab.ypp"
    {
                  end_smodel(NULL);
                  rem_indexing(NULL);
               }
    break;

  case 16:
#line 231 "sml.tab.ypp"
    {
                  end_smodel((char*)(yyvsp[(3) - (4)].string));
                  rem_indexing(NULL);
               }
    break;

  case 18:
#line 245 "sml.tab.ypp"
    {
               stages_glo = NULL;
            }
    break;

  case 19:
#line 253 "sml.tab.ypp"
    {              
                     if (!is_stoch_model) { 
                        cerr << "Syntax Error: keyword 'stages' can only be "
                        "used in stochastic blocks" << endl;
                        exit(1);
                     }
                     stages_glo = (yyvsp[(2) - (3)].opPtr);
                  }
    break;

  case 29:
#line 286 "sml.tab.ypp"
    {(yyval.optype)=LE;}
    break;

  case 30:
#line 287 "sml.tab.ypp"
    {(yyval.optype)=GE;}
    break;

  case 31:
#line 288 "sml.tab.ypp"
    {(yyval.optype)=LT;}
    break;

  case 32:
#line 289 "sml.tab.ypp"
    {(yyval.optype)=GT;}
    break;

  case 33:
#line 290 "sml.tab.ypp"
    {(yyval.optype)=ASSIGN;}
    break;

  case 34:
#line 291 "sml.tab.ypp"
    {(yyval.optype)=EQ;}
    break;

  case 35:
#line 292 "sml.tab.ypp"
    {(yyval.optype)=NE;}
    break;

  case 36:
#line 295 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(ASSIGN, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 37:
#line 296 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(LT, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 38:
#line 297 "sml.tab.ypp"
    { 
               OpNode *temp = new OpNode(LT, (yyvsp[(1) - (5)].opPtr), (yyvsp[(3) - (5)].opPtr));
               (yyval.opPtr) = new OpNode(LT, temp, (yyvsp[(5) - (5)].opPtr)); 
            }
    break;

  case 39:
#line 301 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(LE, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 40:
#line 302 "sml.tab.ypp"
    { 
               OpNode *temp = new OpNode(LE, (yyvsp[(1) - (5)].opPtr), (yyvsp[(3) - (5)].opPtr));
               (yyval.opPtr) = new OpNode(LE, temp, (yyvsp[(5) - (5)].opPtr)); 
            }
    break;

  case 41:
#line 306 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(GT, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 42:
#line 307 "sml.tab.ypp"
    { 
               OpNode *temp = new OpNode(GT, (yyvsp[(1) - (5)].opPtr), (yyvsp[(3) - (5)].opPtr));
               (yyval.opPtr) = new OpNode(GT, temp, (yyvsp[(5) - (5)].opPtr)); 
            }
    break;

  case 43:
#line 311 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(GE, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 44:
#line 312 "sml.tab.ypp"
    { 
               OpNode *temp = new OpNode(GE, (yyvsp[(1) - (5)].opPtr), (yyvsp[(3) - (5)].opPtr));
               (yyval.opPtr) = new OpNode(GE, temp, (yyvsp[(5) - (5)].opPtr)); 
            }
    break;

  case 45:
#line 316 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(EQ, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 46:
#line 319 "sml.tab.ypp"
    {
            if ((yyvsp[(3) - (3)].opPtrIx)) add_indexing((yyvsp[(3) - (3)].opPtrIx));
         }
    break;

  case 47:
#line 321 "sml.tab.ypp"
    {
            ModelComp *newmc;
            if (is_stoch_model){
               newmc = new StochModelComp((yyvsp[(2) - (7)].string), TMAX, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
               addStochInfo(newmc, (yyvsp[(5) - (7)].opPtr));
            }else{
               newmc = new ModelComp((yyvsp[(2) - (7)].string), TMAX, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
            }

            current_model->addComp(newmc);
            current_model->symbol_table.defineSymbol(SymbolTable::ST_OBJ, (yyvsp[(2) - (7)].string),
               newmc);
            if ((yyvsp[(3) - (7)].opPtrIx)) rem_indexing((yyvsp[(3) - (7)].opPtrIx));
            (yyval.string)=(yyvsp[(2) - (7)].string);
         }
    break;

  case 48:
#line 336 "sml.tab.ypp"
    {
            if ((yyvsp[(3) - (3)].opPtrIx)) add_indexing((yyvsp[(3) - (3)].opPtrIx));
         }
    break;

  case 49:
#line 338 "sml.tab.ypp"
    {
            ModelComp *newmc;
            if (is_stoch_model){
               newmc = new StochModelComp((yyvsp[(2) - (7)].string), TMIN, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
               addStochInfo(newmc, (yyvsp[(5) - (7)].opPtr));
            }else{
               newmc = new ModelComp((yyvsp[(2) - (7)].string), TMIN, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
            }

            current_model->addComp(newmc);
            current_model->symbol_table.defineSymbol(SymbolTable::ST_OBJ, (yyvsp[(2) - (7)].string),
               newmc); 
            if ((yyvsp[(3) - (7)].opPtrIx)) rem_indexing((yyvsp[(3) - (7)].opPtrIx));
            (yyval.string)=(yyvsp[(2) - (7)].string);
         }
    break;

  case 50:
#line 353 "sml.tab.ypp"
    {
            if ((yyvsp[(3) - (3)].opPtrIx)) add_indexing((yyvsp[(3) - (3)].opPtrIx));
         }
    break;

  case 51:
#line 355 "sml.tab.ypp"
    { 
            ModelComp *newmc;
            if (is_stoch_model){
               newmc = new StochModelComp((yyvsp[(2) - (7)].string), TCON, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
               addStochInfo(newmc, (yyvsp[(5) - (7)].opPtr));
            }else{
               newmc = new ModelComp((yyvsp[(2) - (7)].string), TCON, (yyvsp[(3) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
            }
            current_model->addComp(newmc);
            current_model->symbol_table.defineSymbol(SymbolTable::ST_CONS, (yyvsp[(2) - (7)].string),
               newmc); 
            if ((yyvsp[(3) - (7)].opPtrIx)) rem_indexing((yyvsp[(3) - (7)].opPtrIx));
            (yyval.string)=(yyvsp[(2) - (7)].string);
         }
    break;

  case 52:
#line 371 "sml.tab.ypp"
    {(yyval.opPtrIx)=NULL;}
    break;

  case 53:
#line 372 "sml.tab.ypp"
    {(yyval.opPtrIx)=(yyvsp[(1) - (1)].opPtrIx);}
    break;

  case 54:
#line 375 "sml.tab.ypp"
    { 
                     SyntaxNodeIx *tmp = 
                        new SyntaxNodeIx(new SyntaxNode(LBRACE, (yyvsp[(2) - (2)].opPtr)));
                     add_indexing(tmp);
                     (yyval.opPtr) = (yyvsp[(2) - (2)].opPtr);
                  }
    break;

  case 55:
#line 383 "sml.tab.ypp"
    {
               rem_indexing(NULL);
               SyntaxNode *tmp = 
                  new SyntaxNode(LBRACE, new SyntaxNode(COLON, (yyvsp[(1) - (4)].opPtr), (yyvsp[(3) - (4)].opPtr)));
               (yyval.opPtrIx) = new SyntaxNodeIx(tmp);
            }
    break;

  case 56:
#line 389 "sml.tab.ypp"
    {
               rem_indexing(NULL);
               SyntaxNodeIx *tmp = new SyntaxNodeIx(new SyntaxNode(LBRACE, (yyvsp[(1) - (2)].opPtr)));
               (yyval.opPtrIx)=tmp;
            }
    break;

  case 58:
#line 400 "sml.tab.ypp"
    {
                  if ((yyvsp[(1) - (3)].opPtr)->getOpCode() == COMMA)
                     (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr)->push_back((yyvsp[(3) - (3)].opPtr));
                  else
                     (yyval.opPtr) = new ListNode(COMMA, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 60:
#line 409 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new OpNode(IN, new IDNode((yyvsp[(1) - (3)].string)), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 61:
#line 412 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new OpNode(IN, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 62:
#line 417 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(LBRACKET, (yyvsp[(1) - (2)].opPtr)); }
    break;

  case 63:
#line 420 "sml.tab.ypp"
    { 
                     (yyval.opPtr) = new ListNode(COMMA, new IDNode((yyvsp[(2) - (2)].string)));
                  }
    break;

  case 64:
#line 423 "sml.tab.ypp"
    {
                     (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr)->push_back(new IDNode((yyvsp[(3) - (3)].string)));
                  }
    break;

  case 66:
#line 429 "sml.tab.ypp"
    {
            (yyval.opPtr) = new SyntaxNode(LBRACKET, (yyvsp[(2) - (3)].opPtr)); }
    break;

  case 68:
#line 434 "sml.tab.ypp"
    {
            if ((yyvsp[(3) - (3)].opPtrIx)) add_indexing((yyvsp[(3) - (3)].opPtrIx));
         }
    break;

  case 69:
#line 436 "sml.tab.ypp"
    {
            ModelComp *newmc;
            if (is_stoch_model){
               newmc = new StochModelComp((yyvsp[(2) - (6)].string), TSET, (yyvsp[(3) - (6)].opPtrIx), (yyvsp[(6) - (6)].opPtr));
               addStochInfo(newmc, (yyvsp[(5) - (6)].opPtr));
            }else{
               newmc = new ModelComp((yyvsp[(2) - (6)].string), TSET, (yyvsp[(3) - (6)].opPtrIx), (yyvsp[(6) - (6)].opPtr));
            }
            current_model->addComp(newmc);
            current_model->symbol_table.defineSymbol(SymbolTable::ST_SET, (yyvsp[(2) - (6)].string),
               newmc); 
            if (GlobalVariables::logParseModel) 
               cout << "$6 = " << (yyvsp[(6) - (6)].opPtr) << "\n";
            if (GlobalVariables::logParseModel)
               cout << "$5 = " << (yyvsp[(5) - (6)].opPtr) << "\n";
            if ((yyvsp[(3) - (6)].opPtrIx)) rem_indexing((yyvsp[(3) - (6)].opPtrIx));
            //$$=$2; 
         }
    break;

  case 70:
#line 457 "sml.tab.ypp"
    {
               if ((yyvsp[(4) - (4)].opPtrIx)) add_indexing((yyvsp[(4) - (4)].opPtrIx));
            }
    break;

  case 71:
#line 459 "sml.tab.ypp"
    {
               ModelComp *newmc;
               if (is_stoch_model){
                  StochModelComp *newmcs = new StochModelComp((yyvsp[(3) - (7)].string),TPARAM,(yyvsp[(4) - (7)].opPtrIx),(yyvsp[(7) - (7)].opPtr));
                  addStochInfo(newmcs, (yyvsp[(6) - (7)].opPtr));
                  newmc = newmcs;
               }else{
                  newmc = new ModelComp((yyvsp[(3) - (7)].string), TPARAM, (yyvsp[(4) - (7)].opPtrIx), (yyvsp[(7) - (7)].opPtr));
               }
               current_model->addComp(newmc);
               current_model->symbol_table.defineSymbol(
                  SymbolTable::ST_PARAM, (yyvsp[(3) - (7)].string), newmc); 
               if ((yyvsp[(4) - (7)].opPtrIx)) rem_indexing((yyvsp[(4) - (7)].opPtrIx));
               //$$=$2; 
            }
    break;

  case 72:
#line 476 "sml.tab.ypp"
    {
            if ((yyvsp[(3) - (3)].opPtrIx)) add_indexing((yyvsp[(3) - (3)].opPtrIx));
         }
    break;

  case 73:
#line 478 "sml.tab.ypp"
    {
            ModelComp *newmc;
            if (is_stoch_model){
               newmc = new StochModelComp((yyvsp[(2) - (6)].string), TVAR, (yyvsp[(3) - (6)].opPtrIx), (yyvsp[(6) - (6)].opPtr));
               addStochInfo(newmc, (yyvsp[(5) - (6)].opPtr));
            }else{
               newmc = new ModelComp((yyvsp[(2) - (6)].string), TVAR, (yyvsp[(3) - (6)].opPtrIx), (yyvsp[(6) - (6)].opPtr));
            }

            current_model->addComp(newmc);
            current_model->symbol_table.defineSymbol(SymbolTable::ST_VAR, (yyvsp[(2) - (6)].string),
               newmc); 
            if ((yyvsp[(3) - (6)].opPtrIx)) rem_indexing((yyvsp[(3) - (6)].opPtrIx));
            //$$=$2; 
         }
    break;

  case 75:
#line 504 "sml.tab.ypp"
    {(yyval.opPtr)=NULL;}
    break;

  case 76:
#line 505 "sml.tab.ypp"
    {(yyval.opPtr)=(yyvsp[(1) - (1)].opPtr);}
    break;

  case 77:
#line 508 "sml.tab.ypp"
    {(yyval.opPtr)=(yyvsp[(1) - (1)].opPtr);}
    break;

  case 78:
#line 509 "sml.tab.ypp"
    {
                  if ((yyvsp[(2) - (2)].opPtr)==NULL){
                     (yyval.opPtr) = (yyvsp[(1) - (2)].opPtr);
                  }else{
                     if ((yyvsp[(1) - (2)].opPtr)==NULL){
                        (yyval.opPtr) = new ListNode(' ', (yyvsp[(2) - (2)].opPtr));
                     }else{
                        if ((yyvsp[(1) - (2)].opPtr)->getOpCode() == ' ')
                           (yyval.opPtr) = (yyvsp[(1) - (2)].opPtr)->push_back((yyvsp[(2) - (2)].opPtr));
                        else
                           (yyval.opPtr) = new ListNode(' ', (yyvsp[(1) - (2)].opPtr), (yyvsp[(2) - (2)].opPtr));
                     }
                  }
               }
    break;

  case 79:
#line 525 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new SyntaxNode(DIMEN, new ValueNode<long>(*(yyvsp[(2) - (2)].ival)));
               }
    break;

  case 80:
#line 528 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new SyntaxNode(WITHIN, (yyvsp[(2) - (2)].opPtr));
               }
    break;

  case 81:
#line 531 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new OpNode(DEFINED, (yyvsp[(2) - (2)].opPtr));
               }
    break;

  case 82:
#line 534 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new SyntaxNode(DEFAULT, (yyvsp[(2) - (2)].opPtr));
               }
    break;

  case 83:
#line 537 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(ORDERED); }
    break;

  case 84:
#line 540 "sml.tab.ypp"
    { (yyval.opPtr) = NULL; }
    break;

  case 85:
#line 541 "sml.tab.ypp"
    { (yyval.opPtr) = (yyvsp[(1) - (1)].opPtr); }
    break;

  case 86:
#line 544 "sml.tab.ypp"
    {
                     if ((yyvsp[(1) - (1)].opPtr)==NULL) {(yyval.opPtr) = NULL;}
                     else{(yyval.opPtr) = new ListNode(COMMA, (yyvsp[(1) - (1)].opPtr));}
                  }
    break;

  case 87:
#line 548 "sml.tab.ypp"
    {
                     if ((yyvsp[(3) - (3)].opPtr)==NULL){
                        (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr);
                     }else{
                        if ((yyvsp[(1) - (3)].opPtr)==NULL){
                           (yyval.opPtr) = new ListNode(COMMA, (yyvsp[(3) - (3)].opPtr));
                        }else{
                           (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr)->push_back((yyvsp[(3) - (3)].opPtr));
                        }
                     }
                  }
    break;

  case 88:
#line 561 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(BINARY);}
    break;

  case 89:
#line 562 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(INTEGER);}
    break;

  case 90:
#line 563 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(SYMBOLIC);}
    break;

  case 91:
#line 564 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode((yyvsp[(1) - (2)].optype), (yyvsp[(2) - (2)].opPtr));}
    break;

  case 92:
#line 565 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(IN, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 93:
#line 566 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(ASSIGN, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 94:
#line 567 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(DEFAULT, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 95:
#line 568 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(DEFINED, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 96:
#line 571 "sml.tab.ypp"
    {(yyval.opPtr) = NULL;}
    break;

  case 97:
#line 572 "sml.tab.ypp"
    {
                  // check that this is in a stochastic model
                  if (!is_stoch_model){ 
                     cerr << "Syntax Error: keyword 'DETERMINISTIC' can only"
                        "be used in stochastic blocks\n";
                     exit(1);
                  }
                  (yyval.opPtr) = new SyntaxNode(DETERMINISTIC);
               }
    break;

  case 98:
#line 581 "sml.tab.ypp"
    {
                  // check that this is in a stochastic model
                  if (!is_stoch_model){ 
                     cerr << "Syntax Error: keyword 'STAGES' can only be used"
                        "in stochastic blocks\n";
                     exit(1);
                  }
                  (yyval.opPtr) = (yyvsp[(2) - (2)].opPtr);
               }
    break;

  case 99:
#line 592 "sml.tab.ypp"
    {(yyval.opPtr)=NULL;}
    break;

  case 100:
#line 593 "sml.tab.ypp"
    {(yyval.opPtr) = (yyvsp[(1) - (1)].opPtr);}
    break;

  case 101:
#line 597 "sml.tab.ypp"
    {
                  (yyval.opPtr) = addItemToListOrCreate(COMMA, NULL, (yyvsp[(1) - (1)].opPtr));}
    break;

  case 102:
#line 599 "sml.tab.ypp"
    {
                  (yyval.opPtr) = addItemToListOrCreate(COMMA, (ListNode*)(yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 103:
#line 605 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(BINARY);}
    break;

  case 104:
#line 606 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(INTEGER);}
    break;

  case 105:
#line 607 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(SYMBOLIC);}
    break;

  case 106:
#line 608 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(LE, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 107:
#line 609 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(GE, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 108:
#line 610 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(DEFINED, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 109:
#line 611 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(ASSIGN, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 110:
#line 612 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(DEFAULT, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 111:
#line 616 "sml.tab.ypp"
    {(yyval.opPtr) = new OpNode(IN, (yyvsp[(2) - (2)].opPtr));}
    break;

  case 112:
#line 617 "sml.tab.ypp"
    {(yyval.opPtr) = new SyntaxNode(SUFFIX, (yyvsp[(3) - (3)].opPtr));}
    break;

  case 113:
#line 621 "sml.tab.ypp"
    {(yyval.opPtr) = new ListSet((yyvsp[(2) - (3)].opPtr));}
    break;

  case 115:
#line 623 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new CompositeSet((yyvsp[(2) - (3)].optype), (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 116:
#line 626 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new CompositeSet((yyvsp[(2) - (3)].optype), (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 117:
#line 629 "sml.tab.ypp"
    {
                  cerr << "FIXME: ubsetop indexing setexpression\n";
                  exit(2);
               }
    break;

  case 118:
#line 633 "sml.tab.ypp"
    {
                  (yyval.opPtr) = new SimpleSet((yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr));
               }
    break;

  case 119:
#line 636 "sml.tab.ypp"
    {
                  cerr << "FIXME: SETOF\n";
                  exit(2);
               }
    break;

  case 120:
#line 641 "sml.tab.ypp"
    {
                  (yyval.opPtr) = (yyvsp[(2) - (3)].opPtr);
               }
    break;

  case 121:
#line 649 "sml.tab.ypp"
    { (yyval.optype) = DIFF; }
    break;

  case 122:
#line 650 "sml.tab.ypp"
    { (yyval.optype) = SYMDIFF; }
    break;

  case 123:
#line 651 "sml.tab.ypp"
    { (yyval.optype) = CROSS; }
    break;

  case 124:
#line 654 "sml.tab.ypp"
    { (yyval.optype) = UNION; }
    break;

  case 125:
#line 655 "sml.tab.ypp"
    { (yyval.optype) = INTER; }
    break;

  case 126:
#line 664 "sml.tab.ypp"
    {
               /* this is a simple identifier in global context */
               (yyval.opPtr) = find_var_ref_in_context(current_model, (yyvsp[(1) - (1)].opPtr));
            }
    break;

  case 127:
#line 668 "sml.tab.ypp"
    {
               /* identifier sets the context for the iditem:
                  The result of this is either a context setting or a
                  complete description of a variable.

                  Can implement this by simply adding context and argument lists
                
                  A variable reference should be represented internally as a
                  pointer to the referenced object in the model list
                  and a pointer to the list of arguments

                  This can be represented as a SyntaxNode with
                  ->opCode = ID;
                  ->nval = # of arguments
                  ->values[0] = pointer to entity in model list
                  ->values[1 - n] = arguments 

                  CONTEXT: iditems can be interpreted in two different contexts
                  1) is the global context (i.e. referring from the model
                     part that is currently defined)
                  2) is a local context that can be set by preceeding bits
                     of a dot'd expression. If a dot'd expression is 
                     parsed the flag 'local_context' should be set and 
                     'local_context' should be set to the current local
                     context. 
                  a context is expressed as a pointer to a model entity
               */
       
               /* identifier sets the context for the idem */
               /* this only works if the identifier is actually a reference 
                  to a submodel */

               if ((yyvsp[(1) - (3)].opPtr)->getOpCode() != IDREFM) {
                  cerr << "Attempting to use dot specifier for something "
                     "not an object:\n " << *((yyvsp[(1) - (3)].opPtr)) << "\n";
                  exit(1);
               }
               local_context = (((SyntaxNodeIDREF*)(yyvsp[(1) - (3)].opPtr))->ref)->other;

               if (GlobalVariables::logParseModel) {
                  cout << "Trying to merge \n identifier "<< *((yyvsp[(1) - (3)].opPtr));
                  cout << " and iditem " << *((yyvsp[(3) - (3)].opPtr)) << "\n";
               }
               (yyval.opPtr) = find_var_ref_in_context(local_context, (yyvsp[(3) - (3)].opPtr));

               /* merge argument lists */
               (yyval.opPtr)->merge((yyvsp[(1) - (3)].opPtr));
            }
    break;

  case 128:
#line 721 "sml.tab.ypp"
    {                               /* simple identifier */
            (yyval.opPtr)=new IDNode((yyvsp[(1) - (1)].string));
         }
    break;

  case 129:
#line 724 "sml.tab.ypp"
    { /* subscripted id'fier */
            if (GlobalVariables::logParseModel) 
               cout << print_SyntaxNodesymb((yyvsp[(3) - (4)].opPtr)) << "\n";
            (yyval.opPtr) = new SyntaxNode(LSBRACKET, new IDNode((yyvsp[(1) - (4)].string)), (yyvsp[(3) - (4)].opPtr));
            //printf("%s\n", print_SyntaxNodesymb($$));
         }
    break;

  case 130:
#line 730 "sml.tab.ypp"
    {
            // This is of the type xh(-1,i) which is xh[i] at a previous stage
            // the ancestor information is conveyed by an ID SyntaxNode with 
            // two arguments, where the second argument is the ancestor
            // $3 is (long int*) => change its sign
            *(yyvsp[(3) - (6)].ival) = -(*(yyvsp[(3) - (6)].ival));
            SyntaxNode *nd = new IDNode((yyvsp[(1) - (6)].string), *(yyvsp[(3) - (6)].ival));
            (yyval.opPtr) = new SyntaxNode(LSBRACKET, nd, (yyvsp[(5) - (6)].opPtr));
         }
    break;

  case 131:
#line 739 "sml.tab.ypp"
    {
            // the same as above, just different syntax "ancestor(1).xh[i]"
            // => need to change the ID node in iditem into a binary node
            //    with INT_VAL ($3) as the second argument
   
            // iditem is either an ID or a LSBRACKET node
            SyntaxNode *node = (yyvsp[(6) - (6)].opPtr);
            IDNode *idnode;
            if (node->getOpCode() == LSBRACKET)
               idnode = (IDNode*) node->front();
            else
               idnode = (IDNode *) node;
            assert(idnode->getOpCode() == ID);
            assert(idnode->getStochParent() == 0);
            idnode->setStochParent(*(yyvsp[(3) - (6)].ival));
            (yyval.opPtr) = (yyvsp[(6) - (6)].opPtr);
         }
    break;

  case 132:
#line 764 "sml.tab.ypp"
    {
               (yyval.opPtr) = new ListNode(COMMA, (yyvsp[(1) - (1)].opPtr));
            }
    break;

  case 133:
#line 767 "sml.tab.ypp"
    {      /* add item to list */
               /* epxr_list could be a simple node or a comma separated
                  list (CSL) already 
                  - if it is a simple node, need to start a comma separated list
                  - if it is a CSL need to add an item to it
               */
               //printf("join nodes >%s< >%s<\n",
               //print_SyntaxNode($1),print_SyntaxNode($3));
               assert((yyvsp[(1) - (3)].opPtr)->getOpCode() == COMMA);
               (yyval.opPtr) = (yyvsp[(1) - (3)].opPtr)->push_back((yyvsp[(3) - (3)].opPtr));
            }
    break;

  case 135:
#line 781 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(LBRACKET, (yyvsp[(2) - (3)].opPtr)); }
    break;

  case 136:
#line 782 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode('+', (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 137:
#line 783 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode('-', (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 138:
#line 784 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode('-', (yyvsp[(2) - (2)].opPtr)); }
    break;

  case 139:
#line 785 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode('*', (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 140:
#line 786 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode('/', (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 141:
#line 787 "sml.tab.ypp"
    { (yyval.opPtr) = new OpNode(POWER, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 142:
#line 788 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(ELLIPSE, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 143:
#line 789 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(LOGICAL_OR, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 144:
#line 790 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(LOGICAL_AND, (yyvsp[(1) - (3)].opPtr), (yyvsp[(3) - (3)].opPtr)); }
    break;

  case 145:
#line 791 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode('!', (yyvsp[(2) - (2)].opPtr)); }
    break;

  case 146:
#line 792 "sml.tab.ypp"
    {add_indexing((yyvsp[(2) - (2)].opPtrIx));}
    break;

  case 147:
#line 792 "sml.tab.ypp"
    { 
         /* reduction operator: do we need to keep track of the ID of the
            dummy variable(s)? */
         (yyval.opPtr) = new SyntaxNode((int)(yyvsp[(1) - (4)].optype), (yyvsp[(2) - (4)].opPtrIx), (yyvsp[(4) - (4)].opPtr));
         rem_indexing((yyvsp[(2) - (4)].opPtrIx));
      }
    break;

  case 148:
#line 798 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(IF, (yyvsp[(2) - (4)].opPtr), (yyvsp[(4) - (4)].opPtr)); }
    break;

  case 149:
#line 799 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(IF, (yyvsp[(2) - (6)].opPtr), (yyvsp[(4) - (6)].opPtr), (yyvsp[(6) - (6)].opPtr)); }
    break;

  case 150:
#line 800 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(FIRST, (yyvsp[(3) - (4)].opPtr)); }
    break;

  case 151:
#line 801 "sml.tab.ypp"
    { (yyval.opPtr) = new SyntaxNode(LAST, (yyvsp[(3) - (4)].opPtr)); }
    break;

  case 152:
#line 802 "sml.tab.ypp"
    { 
         //$$ = new SyntaxNode(EXPECTATION, new SyntaxNode(LBRACKET, $3));}
         (yyval.opPtr) = new SyntaxNode(EXPECTATION, (yyvsp[(3) - (4)].opPtr));
      }
    break;

  case 153:
#line 806 "sml.tab.ypp"
    { /* function definition */
         (yyval.opPtr) = new SyntaxNode((yyvsp[(1) - (4)].optype), new SyntaxNode(LBRACKET, (yyvsp[(3) - (4)].opPtr)));
      }
    break;

  case 154:
#line 811 "sml.tab.ypp"
    { (yyval.optype)=ORD; }
    break;

  case 155:
#line 812 "sml.tab.ypp"
    { (yyval.optype)=CARD; }
    break;

  case 156:
#line 815 "sml.tab.ypp"
    { (yyval.optype)=SUM; }
    break;

  case 157:
#line 816 "sml.tab.ypp"
    { (yyval.optype)=MAX; }
    break;

  case 158:
#line 817 "sml.tab.ypp"
    { (yyval.optype)=MIN; }
    break;

  case 159:
#line 818 "sml.tab.ypp"
    { (yyval.optype)=PROD; }
    break;

  case 160:
#line 821 "sml.tab.ypp"
    {
            (yyval.opPtr)=new ValueNode<long>(*(yyvsp[(1) - (1)].ival));
         }
    break;

  case 161:
#line 824 "sml.tab.ypp"
    { 
            (yyval.opPtr)=new ValueNode<double>(*(yyvsp[(1) - (1)].fval));
         }
    break;

  case 162:
#line 827 "sml.tab.ypp"
    { 
            (yyval.opPtr)=new SyntaxNode(0, (yyvsp[(1) - (1)].opPtr));
         }
    break;

  case 163:
#line 830 "sml.tab.ypp"
    { 
            (yyval.opPtr) = new SyntaxNode(INFINITY);
         }
    break;


/* Line 1267 of yacc.c.  */
#line 3064 "sml.tab.cpp"
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


#line 865 "sml.tab.ypp"


void yyerror(const char *s) {
   cerr << "MODEL: " << s << " on line " << yylineno << "\n";
   //fprintf(stderr, "%s\n", s);
}


/* ----------------------------------------------------------------------------
yywrap
---------------------------------------------------------------------------- */
/* not sure if this is correct, found this somewhere on the internet
   should open the data file and somehow tell the parser to carry on 
   reading in "data mode"

 */
int yywrap(void) {
   return 1;
}

/* ----------------------------------------------------------------------------
begin_model
---------------------------------------------------------------------------- */
void
begin_model(char *name, SyntaxNode *indexing) {
  AmplModel *new_mod = new AmplModel(name);
  ModelComp *newmc;
  
  cout << "Start Model: " << name << "\n";

  /* FIXME: include attrib in definition */
  newmc = new ModelComp(name, TMODEL, indexing, NULL);
  newmc->other = new_mod;

  new_mod->node = newmc;            /* add pointer-to-node to the model */
  current_model->addComp(newmc);
 
  new_mod->setGlobalName();    
  /* and change current model */
  current_model = new_mod;
}

/* ----------------------------------------------------------------------------
begin_smodel
---------------------------------------------------------------------------- */
void
begin_smodel(char *name, SyntaxNode *indexing, SyntaxNode *stochsets) {
   StochModel *new_mod;
   ModelComp *newmc;
  
   ListNode *stochsetsl = static_cast<ListNode*>(stochsets);

   if (!stochsetsl || stochsetsl->nchild()!=4){
      cerr << "Syntax error in Stochastic Block definition: \n";
      cerr << " 'USING' needs 4 parameters \n";
      exit(1);
   }

   if (GlobalVariables::logParseModel)
      cout << "Start Stochastic Model: " << name << "\n";

   ListNode::iterator i = stochsetsl->begin();
   SyntaxNode *nodes = *i;
   SyntaxNode *anc = *(++i);
   SyntaxNode *prob = *(++i);
   SyntaxNode *stages = *(++i);
   /*cout << "BEG SMODEL " << name << endl << "   nodes = " << nodes << endl;
   cout << "   anc = " << anc << endl << "   prob = " << prob << endl;
   cout << "   stages = " << stages << endl;*/
   new_mod = new StochModel(stages, nodes, anc, prob, current_model);
   new_mod->name = name;
  
   /* Fixme: include attrib in definition */
   newmc = new ModelComp(name, TMODEL, indexing, NULL);
   newmc->other = new_mod;

   new_mod->node = newmc;            /* add pointer-to-node to the model */
   current_model->addComp(newmc);
 
   new_mod->setGlobalName();    
   /* and change current model */
   current_model = new_mod;
   is_stoch_model = true;
}

/* ----------------------------------------------------------------------------
end_model
---------------------------------------------------------------------------- */
void
end_model(char *name) {
  // Check end block name matches block name
  if (name && string(name) != current_model->name) {
    cerr << "end block '" << name << "' encountered in block '" << 
      current_model->name << "'" << endl;
    exit(1);
  }

  current_model = current_model->parent;
}

/* ----------------------------------------------------------------------------
end_smodel
---------------------------------------------------------------------------- */
void
end_smodel(char *name){
  // current_model is a StochModel -> convert this into a tree a FlatModels

  // Check end block name matches block name
  if (name && name != current_model->name) {
    cerr << "end stochastic block '" << name << "' encountered in stochastic "
      "block '" << current_model->name << "'" << endl;
    exit(1);
  }

  // this is the ModelComp pointing to the StochModel
  ModelComp *mc = current_model->node; 
  
  // point that to the expanded flat model tree
  mc->other = current_model->expandToFlatModel();

  // and change the name of the ModelComp of this model to the name of the 
  // new (AmplModel) model. 
  // (this is a concatenation of the StochModel name and the name of the 
  // first stage)
  mc->id = mc->other->name;

  // and go back to the parent 
  current_model = current_model->parent;
  is_stoch_model = false;
}

/* ------------------------------------------------------------------------
add_indexing/rem_indexing
-------------------------------------------------------------------------- */
void 
add_indexing(SyntaxNodeIx *indexing){
   list_of_indexing[n_indexing] = indexing;
   if (GlobalVariables::logParseModel){
      cout << "add indexing expression to list: " << *indexing << "\n";
      cout << "Symbolic indexing: " << print_SyntaxNodesymb(indexing) << "\n";
      cout << "length of indexing now: " << n_indexing+1 << "\n";
   }
   n_indexing++;
}

void
rem_indexing(SyntaxNodeIx *indexing){
   if (indexing){
      assert(indexing==list_of_indexing[n_indexing-1]);
      if (GlobalVariables::logParseModel) 
         cout << "rem indexing expression to list: " << *indexing << "\n";
   }
   if (GlobalVariables::logParseModel) 
      cout << "length of indexing now: " << n_indexing-1 << "\n";
   n_indexing--;
}

/* ---------------------------------------------------------------------------
Stochastic model helper functions
---------------------------------------------------------------------------- */
void 
addStochInfo(ModelComp *newmcs, SyntaxNode *stochopt) {
   if(stochopt) {
     bool isDet = stochopt->getOpCode() == DETERMINISTIC;
     newmcs->setDeterministic(isDet || is_deterministic_glo);
     newmcs->setStageSet(isDet ? stages_glo : stochopt);
   } else {
      newmcs->setDeterministic(is_deterministic_glo);
      newmcs->setStageSet(stages_glo);
   }
}

int parse_model(const string& filename) {

   yyin = fopen(filename.c_str(), "r");
   if (!yyin) {
     cerr << "ERROR: Cannot open file '" << filename << "'.\n";
     return 1;
   }

   // check that we can access the datafile, otherwise we get an ugly message
   // from amplsolver in case the file cannot be accessed
   FILE *fin = fopen(GlobalVariables::datafilename.c_str(), "r");
   if (!fin) {
     cerr << "Cannot open file: '" << GlobalVariables::datafilename << "'.\n";
     return 1;
   }
   fclose(fin);

   // Change directory to make things work
   int errcode = chdir("tmp");
   if (errcode){
      cerr << "Could not change working directory to 'tmp/'\n";
      cerr << "Cannot continue\n";
      return 1;
   }


   AmplModel::root = new AmplModel("root");
   current_model = AmplModel::root;
   is_stoch_model = false;
   is_deterministic_glo = false;
   stages_glo = NULL;

   errcode = yyparse();
   if (errcode!=0){
      cerr << "yyparse returns " << errcode << endl;
      return 1;
   } 

   fclose(yyin);

   // Restore original directory
   errcode = chdir("..");
   if (errcode){
      cerr << "Could not change working directory to 'tmp/'\n";
      cerr << "Cannot continue\n";
      return 1;
   }

   return 0;
}

