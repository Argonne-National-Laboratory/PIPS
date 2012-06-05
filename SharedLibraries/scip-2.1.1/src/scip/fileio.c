/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   fileio.c
 * @brief  wrapper functions to map file i/o to standard or zlib file i/o
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdarg.h>

#include "scip/pub_fileio.h"


#define BUFFER_LEN 8192

#ifdef WITH_ZLIB

/* file i/o using zlib */
#include <zlib.h>

struct SCIP_File
{
   gzFile                file;               /**< zlib file stream */
};

SCIP_FILE* SCIPfopen(const char *path, const char *mode)
{
   return (SCIP_FILE*)gzopen(path, mode);
}

SCIP_FILE* SCIPfdopen(int fildes, const char *mode)
{
   return (SCIP_FILE*)gzdopen(fildes, mode);
}

size_t SCIPfread(void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream)
{
   return gzread((gzFile*)stream, ptr, size * nmemb);
}

size_t SCIPfwrite(const void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream)
{
   return gzwrite((gzFile*)stream, ptr, size * nmemb);
}

int SCIPfprintf(SCIP_FILE *stream, const char *format, ...)
{
   char buffer[BUFFER_LEN];
   va_list ap;
   int n;

   va_start(ap, format); /*lint !e826*/
   n = vsnprintf(buffer, BUFFER_LEN, format, ap);
   va_end(ap);
   if( n < 0 || n > BUFFER_LEN)
      buffer[BUFFER_LEN-1] = '\0';
   return gzputs((gzFile*)stream, buffer);
}

int SCIPfputc(int c, SCIP_FILE *stream)
{
   return gzputc((gzFile*)stream, c);
}

int SCIPfputs(const char *s, SCIP_FILE *stream)
{
   return gzputs((gzFile*)stream, s);
}

int SCIPfgetc(SCIP_FILE *stream)
{
   return gzgetc((gzFile*)stream);
}

char* SCIPfgets(char *s, int size, SCIP_FILE *stream)
{
   return gzgets((gzFile*)stream, s, size);
}

int SCIPfflush(SCIP_FILE *stream)
{
   return gzflush((gzFile*)stream, Z_SYNC_FLUSH);
}

int SCIPfseek(SCIP_FILE *stream, long offset, int whence)
{
   return gzseek((gzFile*)stream, offset, whence);
}

void SCIPrewind(SCIP_FILE *stream)
{
   gzrewind((gzFile*)stream);
}

long SCIPftell(SCIP_FILE *stream)
{
   return gztell((gzFile*)stream);
}

int SCIPfeof(SCIP_FILE *stream)
{
   return gzeof((gzFile*)stream);
}

int SCIPfclose(SCIP_FILE *fp)
{
   return gzclose((gzFile*)fp);
}


#else


/* file i/o using standard i/o */

SCIP_FILE* SCIPfopen(const char *path, const char *mode)
{
   return (SCIP_FILE*)fopen(path, mode);
}

SCIP_FILE* SCIPfdopen(int fildes, const char *mode)
{
   return (SCIP_FILE*)fdopen(fildes, mode);
}

size_t SCIPfread(void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream)
{
   return fread(ptr, size, nmemb, (FILE*)stream);
}

size_t SCIPfwrite(const void *ptr, size_t size, size_t nmemb, SCIP_FILE *stream)
{
   return fwrite(ptr, size, nmemb, (FILE*)stream);
}

int SCIPfprintf(SCIP_FILE *stream, const char *format, ...)
{
   va_list ap;
   int retval;

   va_start(ap, format); /*lint !e826*/
   retval = vfprintf((FILE*)stream, format, ap);
   va_end(ap);

   return retval;
}

int SCIPfputc(int c, SCIP_FILE *stream)
{
   return fputc(c, (FILE*)stream);
}

int SCIPfputs(const char *s, SCIP_FILE *stream)
{
   return fputs(s, (FILE*)stream);
}

int SCIPfgetc(SCIP_FILE *stream)
{
   return fgetc((FILE*)stream);
}

char* SCIPfgets(char *s, int size, SCIP_FILE *stream)
{
   return fgets(s, size, (FILE*)stream);
}

int SCIPfflush(SCIP_FILE *stream)
{
   return fflush((FILE*)stream);
}

int SCIPfseek(SCIP_FILE *stream, long offset, int whence)
{
   return fseek((FILE*)stream, offset, whence);
}

void SCIPrewind(SCIP_FILE *stream)
{
   rewind((FILE*)stream);
}

long SCIPftell(SCIP_FILE *stream)
{
   return ftell((FILE*)stream);
}

int SCIPfeof(SCIP_FILE *stream)
{
   return feof((FILE*)stream);
}

int SCIPfclose(SCIP_FILE *fp)
{
   return fclose((FILE*)fp);
}


#endif
