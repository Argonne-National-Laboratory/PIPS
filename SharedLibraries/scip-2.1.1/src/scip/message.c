/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This1 file is part of the program and library             */
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

/**@file   message.c
 * @brief  message output methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifndef NPARASCIP
#include <pthread.h>
#include <time.h>
#endif

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/message.h"
#include "scip/misc.h"


#ifndef va_copy
#define va_copy(dest, src) do { BMScopyMemory(&dest, &src); } while( 0 )
#endif

/** error message print method of default message handler */
static
SCIP_DECL_MESSAGEERROR(messageErrorDefault)
{  /*lint --e{715}*/
   fputs(msg, file);
   fflush(file);
}

/** warning message print method of default message handler */
static
SCIP_DECL_MESSAGEWARNING(messageWarningDefault)
{  /*lint --e{715}*/
   fputs(msg, file);
   fflush(file);
}

/** dialog message print method of default message handler */
static
SCIP_DECL_MESSAGEDIALOG(messageDialogDefault)
{  /*lint --e{715}*/
   fputs(msg, file);
   fflush(file);
}

/** info message print method of default message handler */
static
SCIP_DECL_MESSAGEINFO(messageInfoDefault)
{  /*lint --e{715}*/
   fputs(msg, file);
   fflush(file);
}

/** default message handler that prints messages to stdout or stderr */
static
SCIP_MESSAGEHDLR messagehdlrDefault =
   { messageErrorDefault, messageWarningDefault, messageDialogDefault, messageInfoDefault,
     NULL, NULL, NULL, NULL, NULL, 0, 0, 0, 0 };


/** stores the given message in the buffer and returns the part of the buffer and message that should be printed now */
static
SCIP_Bool bufferMessage(
   char*                 buffer,             /**< message buffer */
   int*                  bufferlen,          /**< pointer to the currently used entries in the message buffer */
   const char*           msg,                /**< message to store in the buffer; NULL to flush the buffer */
   char*                 outmsg,             /**< array to store message that should be printed immediately */
   int*                  outmsgsize          /**< pointer which stored allocated space in outmsg, size should be at
                                              *   least SCIP_MAXSTRLEN, if space would be not enough the needed sized is
                                              *   returned
                                              */
   )
{
   assert(outmsg != NULL);
   assert(outmsgsize != NULL);
   assert(*outmsgsize >= SCIP_MAXSTRLEN);

   *outmsg = '\0';

   /* should the buffer be flushed? */
   if( msg == NULL )
   {
      if( buffer != NULL )
         (void)strncpy(outmsg, buffer, SCIP_MAXSTRLEN);
      (*bufferlen) = 0;
      assert((int)strlen(outmsg) < SCIP_MAXSTRLEN);
      return TRUE;
   }

   /* do we have a buffer? */
   if( buffer == NULL )
   {
      /* check if the message fits into the output message */
      if( (int)strlen(msg) >= *outmsgsize )
      {
         *outmsgsize = (int)strlen(msg) + 1;
         return FALSE;
      }
      /* no buffer exists -> just copy the message to the output */
      (void)strncpy(outmsg, msg, strlen(msg));
      outmsg[strlen(msg)] = '\0';
      assert((int)strlen(outmsg) < *outmsgsize);
   }
   else
   {
      char* outmsgpos;

      assert(bufferlen != NULL);

      /* check if the message fits into the output message */
      if( (int)strlen(msg) + *bufferlen >= *outmsgsize + SCIP_MAXSTRLEN-1 )
      {
         *outmsgsize = (int)strlen(msg) + *bufferlen - SCIP_MAXSTRLEN + 2;
         return FALSE;
      }

      outmsgpos = outmsg;

      while( *msg != '\0' )
      {
         char c;

         assert(*bufferlen < SCIP_MAXSTRLEN-1);
         c = *msg;
         msg++;
         buffer[*bufferlen] = c;
         (*bufferlen)++;

         /* if the buffer is full or we reached a newline, move the buffer to the output message and store the
          * remaining message in the buffer
          */
         if( *bufferlen >= SCIP_MAXSTRLEN-1 || c == '\n' )
         {
            buffer[*bufferlen] = '\0';
            (void)strncpy(outmsgpos, buffer, SCIP_MAXSTRLEN);
            outmsgpos = &(outmsgpos[*bufferlen]);

            /* if rest fits into buffer copy all */
            if( (int)strlen(msg) < SCIP_MAXSTRLEN - 1 )
            {
               (void)strncpy(buffer, msg, SCIP_MAXSTRLEN);
               *bufferlen = (int)strlen(msg);
               assert(*bufferlen < SCIP_MAXSTRLEN-1);
               assert((int)strlen(outmsg) < *outmsgsize);
               break;
            }
            else
            {
               (void)strncpy(buffer, msg, SCIP_MAXSTRLEN - 2);
               *bufferlen = SCIP_MAXSTRLEN - 2;
               msg += (SCIP_MAXSTRLEN - 2);
            }
         }
      }
   }

   return TRUE;
}


/** for parascip we need to make this threadsafe */
#ifdef NPARASCIP


/** static variable that contains the currently installed message handler;
 *  if the handler is set to NULL, messages are suppressed
 */
static SCIP_MESSAGEHDLR* curmessagehdlr = &messagehdlrDefault;


/** prints error message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintError(
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messageerror != NULL )
   {
      char* outmsg;
      int outmsgsize;

      outmsgsize = msglength + 1;

      if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
      {
         return;
      }

      if( !bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize) )
      {
         assert(outmsgsize > msglength + 1);
         if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
         {
#ifndef NDEBUG
            SCIP_Bool ret;

            ret = bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize);
            assert(ret);
#else
            bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize);
#endif
         }
         else
         {
            BMSfreeMemory(&outmsg);
            return;
         }
      }

      if( *outmsg != '\0' )
         curmessagehdlr->messageerror(curmessagehdlr, stderr, outmsg);

      BMSfreeMemory(&outmsg);
   }
}

/** prints warning message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintWarning(
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messagewarning != NULL )
   {
      char* outmsg;
      int outmsgsize;

      outmsgsize = msglength + 1;

      if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
      {
         return;
      }

      if( !bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize) )
      {
         assert(outmsgsize > msglength + 1);
         if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
         {
#ifndef NDEBUG
            SCIP_Bool ret;

            ret = bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
            assert(ret);
#else
            bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
#endif
         }
         else
         {
            BMSfreeMemory(&outmsg);
            return;
         }
      }

      if( *outmsg != '\0' )
         curmessagehdlr->messagewarning(curmessagehdlr, stderr, outmsg);

      BMSfreeMemory(&outmsg);
   }
}

/** prints dialog message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messagedialog != NULL )
   {
      if( file == NULL || file == stdout )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
         {
            return;
         }

         if( !bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
            curmessagehdlr->messagedialog(curmessagehdlr, stdout, outmsg);

         BMSfreeMemory(&outmsg);
      }
      else
      {
         /* file output cannot be buffered because the output file may change */
         if( *msg != '\0' )
            curmessagehdlr->messagedialog(curmessagehdlr, file, msg);
      }
   }
}

/** prints info message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                 /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   if( curmessagehdlr != NULL && curmessagehdlr->messageinfo != NULL )
   {
      if( (file == NULL || file == stdout) && (msg == NULL || strlen(msg) < SCIP_MAXSTRLEN) )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
         {
            return;
         }

         if( !bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
            curmessagehdlr->messageinfo(curmessagehdlr, stdout, outmsg);

         BMSfreeMemory(&outmsg);
      }
      else
      {
         /* file output cannot be buffered because the output file may change or the message is to long */
         if( *msg != '\0' )
            curmessagehdlr->messageinfo(curmessagehdlr, file == NULL ? stdout : file, msg);
      }
   }
}

/** installs the given message handler */
void SCIPmessageSetHandler(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   curmessagehdlr = messagehdlr;
}

/** installs the default message handler that prints messages to stdout or stderr */
void SCIPmessageSetDefaultHandler(
   void
   )
{
   curmessagehdlr = &messagehdlrDefault;
}

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
SCIP_MESSAGEHDLR* SCIPmessageGetHandler(
   void
   )
{
   return curmessagehdlr;
}


#else /* NPARASCIP */

/* mutex to lock all message printings in pthread case */
static pthread_mutex_t  messagemutex = PTHREAD_MUTEX_INITIALIZER;

/* hashmap which identifies for each thread the corresponding messagehandler */
static SCIP_HASHMAP* messagepthreadhashmap = NULL;

/** static array that should contain all messagehandlers in multiple thread case;
 *  if the handler is set to NULL, messages are suppressed;
 *
 *  @note: in this case it is mandatory to call SCIPcreateMesshdlrPThreads
 */
static SCIP_MESSAGEHDLR** curmessagehdlrs = NULL;
static size_t ncurmessagehdlrs = 0;

/** number of threads for which curmessagehdlrs allocates space, i.e., size of curmessagehdlrs array */
static size_t nthreads = 0;

/* standard message handler when not initializing the array above by calling */
static SCIP_MESSAGEHDLR* emergencycurmessagehdlr = &messagehdlrDefault;

/** prints error message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintError(
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   pthread_mutex_lock(&messagemutex);

   {
      SCIP_MESSAGEHDLR* curmessagehdlr;

      if( ncurmessagehdlrs == 0 )
         curmessagehdlr = emergencycurmessagehdlr;
      else
      {
         assert(messagepthreadhashmap != NULL);
         assert(curmessagehdlrs != NULL);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);
         
         curmessagehdlr = curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1];
      }

      if( curmessagehdlr != NULL && curmessagehdlr->messageerror != NULL )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
         {
            return;
         }

         if( !bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(curmessagehdlr->errorbuffer, &curmessagehdlr->errorbufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
            curmessagehdlr->messageerror(curmessagehdlr, stderr, outmsg);

         BMSfreeMemory(&outmsg);
      }
   }

   pthread_mutex_unlock(&messagemutex);   
}

/** prints warning message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintWarning(
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   pthread_mutex_lock(&messagemutex);

   {
      SCIP_MESSAGEHDLR* curmessagehdlr;

      if( ncurmessagehdlrs == 0 )
         curmessagehdlr = emergencycurmessagehdlr;
      else
      {
         assert(messagepthreadhashmap != NULL);
         assert(curmessagehdlrs != NULL);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

         curmessagehdlr = curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1];
      }

      if( curmessagehdlr != NULL && curmessagehdlr->messagewarning != NULL )
      {
         char* outmsg;
         int outmsgsize;

         outmsgsize = msglength + 1;

         if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
         {
            return;
         }

         if( !bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize) )
         {
            assert(outmsgsize > msglength + 1);
            if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
            {
#ifndef NDEBUG
               SCIP_Bool ret;

               ret = bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
               assert(ret);
#else
               bufferMessage(curmessagehdlr->warningbuffer, &curmessagehdlr->warningbufferlen, msg, outmsg, &outmsgsize);
#endif
            }
            else
            {
               BMSfreeMemory(&outmsg);
               return;
            }
         }

         if( *outmsg != '\0' )
            curmessagehdlr->messagewarning(curmessagehdlr, stderr, outmsg);

         BMSfreeMemory(&outmsg);
      }
   }

   pthread_mutex_unlock(&messagemutex);
}

/** prints dialog message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   pthread_mutex_lock(&messagemutex);

   {
      SCIP_MESSAGEHDLR* curmessagehdlr;

      if( ncurmessagehdlrs == 0 )
         curmessagehdlr = emergencycurmessagehdlr;
      else
      {
         assert(messagepthreadhashmap != NULL);
         assert(curmessagehdlrs != NULL);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

         curmessagehdlr = curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1];
      }

      if( curmessagehdlr != NULL && curmessagehdlr->messagedialog != NULL )
      {
         if( file == NULL || file == stdout )
         {
            char* outmsg;
            int outmsgsize;

            outmsgsize = msglength + 1;

            if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
            {
               return;
            }

            if( !bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize) )
            {
               assert(outmsgsize > msglength + 1);
               if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
               {
#ifndef NDEBUG
                  SCIP_Bool ret;

                  ret = bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
                  assert(ret);
#else
                  bufferMessage(curmessagehdlr->dialogbuffer, &curmessagehdlr->dialogbufferlen, msg, outmsg, &outmsgsize);
#endif
               }
               else
               {
                  BMSfreeMemory(&outmsg);
                  return;
               }
            }

            if( *outmsg != '\0' )
               curmessagehdlr->messagedialog(curmessagehdlr, stdout, outmsg);

            BMSfreeMemory(&outmsg);
         }
         else
         {
            /* file output cannot be buffered because the output file may change */
            if( *msg != '\0' )
               curmessagehdlr->messagedialog(curmessagehdlr, file, msg);
         }
      }
   }

   pthread_mutex_unlock(&messagemutex);
}

/** prints info message with the current message handler, or buffers the message if no newline exists */
static
void messagePrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           msg,                /**< message to print; NULL to flush the output buffer */
   int                   msglength           /**< message length if bigger than SCIP_MAXSTRLEN, or SCIP_MAXSTRLEN */
   )
{
   pthread_mutex_lock(&messagemutex);

   {
      SCIP_MESSAGEHDLR* curmessagehdlr;

      if( ncurmessagehdlrs == 0 )
         curmessagehdlr = emergencycurmessagehdlr;
      else
      {
         assert(messagepthreadhashmap != NULL);
         assert(curmessagehdlrs != NULL);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
         assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

         curmessagehdlr = curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1];
      }

      if( curmessagehdlr != NULL && curmessagehdlr->messageinfo != NULL )
      {
         if( (file == NULL || file == stdout) && (msg == NULL || strlen(msg) < SCIP_MAXSTRLEN) )
         {
            char* outmsg;
            int outmsgsize;

            outmsgsize = msglength + 1;

            if( BMSallocMemorySize(&outmsg, outmsgsize) == NULL )
            {
               return;
            }

            if( !bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize) )
            {
               assert(outmsgsize > msglength + 1);
               if( BMSreallocMemorySize(&outmsg, outmsgsize) != NULL )
               {
#ifndef NDEBUG
                  SCIP_Bool ret;

                  ret = bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
                  assert(ret);
#else
                  bufferMessage(curmessagehdlr->infobuffer, &curmessagehdlr->infobufferlen, msg, outmsg, &outmsgsize);
#endif
               }
               else
               {
                  BMSfreeMemory(&outmsg);
                  return;
               }
            }

            if( *outmsg != '\0' )
               curmessagehdlr->messageinfo(curmessagehdlr, stdout, outmsg);

            BMSfreeMemory(&outmsg);
         }
         else
         {
            /* file output cannot be buffered because the output file may change or the message is to long */
            if( *msg != '\0' )
               curmessagehdlr->messageinfo(curmessagehdlr, file == NULL ? stdout : file, msg);
         }
      }
   }

   pthread_mutex_unlock(&messagemutex);
}

/** installs the given message handler */
void SCIPmessageSetHandler(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to install, or NULL to suppress all output */
   )
{
   pthread_mutex_lock(&messagemutex);

   assert(messagepthreadhashmap != NULL);
   assert(curmessagehdlrs != NULL);

   if( SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self()) == NULL )
   {
      ++ncurmessagehdlrs;
      SCIP_CALL_ABORT( SCIPhashmapInsert(messagepthreadhashmap, (void*) pthread_self(), (void*) ncurmessagehdlrs) );
   }   
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

   curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1] = messagehdlr;

   pthread_mutex_unlock(&messagemutex);
}

/** installs the default message handler that prints messages to stdout or stderr */
void SCIPmessageSetDefaultHandler(
   void
   )
{
   pthread_mutex_lock(&messagemutex);

   assert(messagepthreadhashmap != NULL);
   assert(curmessagehdlrs != NULL);

   if( SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self()) == NULL )
   {
      ++ncurmessagehdlrs;
      SCIP_CALL_ABORT( SCIPhashmapInsert(messagepthreadhashmap, (void*) pthread_self(), (void*) ncurmessagehdlrs) );
   }
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);
   
   curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1] = &messagehdlrDefault;

   pthread_mutex_unlock(&messagemutex);
}

/** returns the currently installed message handler, or NULL if messages are currently suppressed */
SCIP_MESSAGEHDLR* SCIPmessageGetHandler(
   void
   )
{
   SCIP_MESSAGEHDLR* curmessagehdlr;

   pthread_mutex_unlock(&messagemutex);

   assert(messagepthreadhashmap != NULL);
   assert(curmessagehdlrs != NULL);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

   curmessagehdlr = curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1];

   pthread_mutex_unlock(&messagemutex);

   return curmessagehdlr;
}

#define HASHMAPSIZE_FACTOR 5

/** allocates memory for all message handlers for number of given threads */
SCIP_RETCODE SCIPmesshdlrCreatePThreads(
   int                   nthreads_           /**< number of threads to allocate memory for */
   )
{
   pthread_mutex_lock(&messagemutex);

   assert(nthreads_ > 0);
   assert(messagepthreadhashmap == NULL);
   assert(curmessagehdlrs == NULL);
   assert(nthreads == 0);
   assert(ncurmessagehdlrs == 0);

   /* check that we only create the hashmap and array once */
   if( messagepthreadhashmap != NULL || curmessagehdlrs != NULL )
   {
      return SCIP_INVALIDCALL;
   }

   nthreads = nthreads_;

   /* create hashmap and message handler array */
   SCIP_CALL( SCIPhashmapCreate(&messagepthreadhashmap, NULL, HASHMAPSIZE_FACTOR * nthreads) );
   SCIP_ALLOC( BMSallocMemoryArray(&curmessagehdlrs, nthreads) );

   pthread_mutex_unlock(&messagemutex);

   return SCIP_OKAY;
}

/** frees memory for all message handlers */
void SCIPmesshdlrFreePThreads(
   void
   )
{
   pthread_mutex_lock(&messagemutex);

   assert(messagepthreadhashmap != NULL);
   assert(curmessagehdlrs != NULL);

   /* free hashmap and message handler array */
   SCIPhashmapFree(&messagepthreadhashmap);
   BMSfreeMemoryArray(&curmessagehdlrs);
   ncurmessagehdlrs = 0;
   nthreads = 0;

   pthread_mutex_unlock(&messagemutex);
}

/** gets maximal number of threads for which message handler space has been allocated
 * returned value is the one given by SCIPmesshdlrCreatePThreads(), or 0 if that hasn't been called
 */
extern
size_t SCIPmessagehdlrGetNThreads(
   void
   )
{
   return nthreads;
}

/** gets number of current thread as assigned by message handler hashmap
 * assigns a new number to the thread if new */
extern
size_t SCIPmessagehdlrGetThreadNum(
   void
   )
{
   if( nthreads == 0 )
      return 0;

   assert(messagepthreadhashmap != NULL);

   if( SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self()) == NULL )
   {
      pthread_mutex_lock(&messagemutex);

      ++ncurmessagehdlrs;
      SCIP_CALL_ABORT( SCIPhashmapInsert(messagepthreadhashmap, (void*) pthread_self(), (void*) ncurmessagehdlrs) );

      assert(ncurmessagehdlrs <= nthreads);
      curmessagehdlrs[((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) - 1] = &messagehdlrDefault;

      pthread_mutex_unlock(&messagemutex);
   }
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) > 0);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= nthreads);
   assert(((size_t) SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self())) <= ncurmessagehdlrs);

   return (size_t)SCIPhashmapGetImage(messagepthreadhashmap, (void*) pthread_self()) - 1;
}

#endif /* NPARASCIP */


/** creates a message handler */
SCIP_RETCODE SCIPmessagehdlrCreate(
   SCIP_MESSAGEHDLR**    messagehdlr,        /**< pointer to store the message handler */
   SCIP_Bool             bufferedoutput,     /**< should the output be buffered up to the next newline? */
   SCIP_DECL_MESSAGEERROR((*messageerror)),  /**< error message print method of message handler */
   SCIP_DECL_MESSAGEWARNING((*messagewarning)),/**< warning message print method of message handler */
   SCIP_DECL_MESSAGEDIALOG((*messagedialog)),/**< dialog message print method of message handler */
   SCIP_DECL_MESSAGEINFO ((*messageinfo)),   /**< info message print method of message handler */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< message handler data */
   )
{
   SCIP_ALLOC( BMSallocMemory(messagehdlr) );
   (*messagehdlr)->messageerror = messageerror;
   (*messagehdlr)->messagewarning = messagewarning;
   (*messagehdlr)->messagedialog = messagedialog;
   (*messagehdlr)->messageinfo = messageinfo;
   (*messagehdlr)->messagehdlrdata = messagehdlrdata;
   (*messagehdlr)->errorbuffer = NULL;
   (*messagehdlr)->warningbuffer = NULL;
   (*messagehdlr)->dialogbuffer = NULL;
   (*messagehdlr)->infobuffer = NULL;
   (*messagehdlr)->errorbufferlen = 0;
   (*messagehdlr)->warningbufferlen = 0;
   (*messagehdlr)->dialogbufferlen = 0;
   (*messagehdlr)->infobufferlen = 0;

   /* allocate buffer for buffered output */
   if( bufferedoutput )
   {
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->errorbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->warningbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->dialogbuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      SCIP_ALLOC( BMSallocMemoryArray(&(*messagehdlr)->infobuffer, SCIP_MAXSTRLEN) ); /*lint !e506*/
      (*messagehdlr)->errorbuffer[0] = '\0';
      (*messagehdlr)->warningbuffer[0] = '\0';
      (*messagehdlr)->dialogbuffer[0] = '\0';
      (*messagehdlr)->infobuffer[0] = '\0';
   }

   return SCIP_OKAY;
}

/** frees message handler */
void SCIPmessagehdlrFree(
   SCIP_MESSAGEHDLR**    messagehdlr         /**< pointer to the message handler */
   )
{
   assert(messagehdlr != NULL);

   if( *messagehdlr != NULL )
   {
      /* flush message buffers */
      messagePrintError(NULL, SCIP_MAXSTRLEN);
      messagePrintWarning(NULL, SCIP_MAXSTRLEN);
      messagePrintDialog(NULL, NULL, SCIP_MAXSTRLEN);
      messagePrintInfo(NULL, NULL, SCIP_MAXSTRLEN);

      BMSfreeMemoryArrayNull(&(*messagehdlr)->errorbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->warningbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->dialogbuffer);
      BMSfreeMemoryArrayNull(&(*messagehdlr)->infobuffer);
      BMSfreeMemory(messagehdlr);
   }
}

/** sets the user data of the message handler */
SCIP_RETCODE SCIPmessagehdlrSetData(
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler; must not be NULL */
   SCIP_MESSAGEHDLRDATA* messagehdlrdata     /**< new message handler data to attach to the handler */
   )
{
   if( messagehdlr == NULL )
   {
      SCIPerrorMessage("cannot set message handler data - no message handler present\n");
      return SCIP_INVALIDDATA;
   }

   messagehdlr->messagehdlrdata = messagehdlrdata;

   return SCIP_OKAY;
}

/** returns the user data of the message handler */
SCIP_MESSAGEHDLRDATA* SCIPmessagehdlrGetData(
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler */
   )
{
   if( messagehdlr != NULL )
      return messagehdlr->messagehdlrdata;
   else
      return NULL;
}

/** prints the header with source file location for an error message */
void SCIPmessagePrintErrorHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[SCIP_MAXSTRLEN];

   /* safe string printing - do not use SCIPsnprintf() since message.c should be independent */
   (void) snprintf(msg, SCIP_MAXSTRLEN, "[%s:%d] ERROR: ", sourcefile, sourceline);
   msg[SCIP_MAXSTRLEN-1] = '\0';
   messagePrintError(msg, SCIP_MAXSTRLEN);
}

/** prints an error message, acting like the printf() command */
void SCIPmessagePrintError(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   char msg[SCIP_MAXSTRLEN];
   va_list ap;
   int n;

   va_start(ap, formatstr); /*lint !e826*/
   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap); /*lint !e718 !e746*/
   va_end(ap);

   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif
      va_list aq;
   
      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
         return;
      
      va_start(aq, formatstr); /*lint !e826*/
#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e718 !e746*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq); /*lint !e718 !e746*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintError(bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }   

   messagePrintError(msg, SCIP_MAXSTRLEN);
}

/** prints the header with source file location for an error message */
void SCIPmessagePrintWarningHeader(
   const char*           sourcefile,         /**< name of the source file that called the function */
   int                   sourceline          /**< line in the source file where the function was called */
   )
{
   char msg[SCIP_MAXSTRLEN];

   /* safe string printing - do not use SCIPsnprintf() since message.c should be independent */
   (void) snprintf(msg, SCIP_MAXSTRLEN, "[%s:%d] Warning: ", sourcefile, sourceline);
   msg[SCIP_MAXSTRLEN-1] = '\0';   
   messagePrintWarning(msg, SCIP_MAXSTRLEN);
}

/** prints a warning message, acting like the printf() command */
void SCIPmessagePrintWarning(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   char msg[SCIP_MAXSTRLEN];
   va_list ap;
   int n;

   va_start(ap, formatstr); /*lint !e826*/
   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   va_end(ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      va_list aq;
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
         return;

      va_start(aq, formatstr); /*lint !e826*/
#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);  /*lint !e718 !e746*/
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);  /*lint !e718 !e746*/
#endif
      assert(m == n);
      va_end(aq);
      messagePrintWarning(bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }   
   messagePrintWarning(msg, SCIP_MAXSTRLEN);
}

/** prints a dialog message that requests user interaction, acting like the printf() command */
void SCIPmessagePrintDialog(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(NULL, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction, acting like the vprintf() command */
void SCIPmessageVPrintDialog(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintDialog(NULL, formatstr, ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the fprintf() command */
void SCIPmessageFPrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintDialog(file, formatstr, ap);
   va_end(ap);
}

/** prints a dialog message that requests user interaction into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintDialog(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;
   
   va_copy(aq, ap);

   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
      assert(m == n);
      va_end(aq);
      messagePrintDialog(file, bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }   
   messagePrintDialog(file, msg, SCIP_MAXSTRLEN);
   va_end(aq);
}

/** prints a message, acting like the printf() command */
void SCIPmessagePrintInfo(
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message, acting like the vprintf() command */
void SCIPmessageVPrintInfo(
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintInfo(NULL, formatstr, ap);
}

/** prints a message into a file, acting like the fprintf() command */
void SCIPmessageFPrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintInfo(file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file, acting like the vfprintf() command */
void SCIPmessageVFPrintInfo(
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   char msg[SCIP_MAXSTRLEN];
   int n;
   va_list aq;

   va_copy(aq, ap);
   
   n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap); 
   if( n < 0 )
      msg[SCIP_MAXSTRLEN-1] = '\0';
   else if( n >= SCIP_MAXSTRLEN )
   {
      char* bigmsg;
#ifndef NDEBUG
      int m;
#endif

      if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
      {
         va_end(aq);
         return;
      }

#ifndef NDEBUG
      m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
      vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
      assert(m == n);
      va_end(aq);
      messagePrintInfo(file, bigmsg, n);
      BMSfreeMemory(&bigmsg);
      return;
   }   
   messagePrintInfo(file, msg, SCIP_MAXSTRLEN);
   va_end(aq);
}

/** prints a message depending on the verbosity level, acting like the printf() command */
void SCIPmessagePrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, NULL, formatstr, ap);
   va_end(ap);
}

/** prints a message depending on the verbosity level, acting like the vprintf() command */
void SCIPmessageVPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, NULL, formatstr, ap);
}

/** prints a message into a file depending on the verbosity level, acting like the fprintf() command */
void SCIPmessageFPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   ...                                       /**< format arguments line in printf() function */
   )
{
   va_list ap;

   va_start(ap, formatstr); /*lint !e826*/
   SCIPmessageVFPrintVerbInfo(verblevel, msgverblevel, file, formatstr, ap);
   va_end(ap);
}

/** prints a message into a file depending on the verbosity level, acting like the vfprintf() command */
void SCIPmessageVFPrintVerbInfo(
   SCIP_VERBLEVEL        verblevel,          /**< current verbosity level */
   SCIP_VERBLEVEL        msgverblevel,       /**< verbosity level of this message */
   FILE*                 file,               /**< file stream to print into, or NULL for stdout */
   const char*           formatstr,          /**< format string like in printf() function */
   va_list               ap                  /**< variable argument list */
   )
{
   assert(msgverblevel > SCIP_VERBLEVEL_NONE);
   assert(msgverblevel <= SCIP_VERBLEVEL_FULL);
   assert(verblevel <= SCIP_VERBLEVEL_FULL);

   if( msgverblevel <= verblevel )
   {
      char msg[SCIP_MAXSTRLEN];
      int n;
      va_list aq;

      va_copy(aq, ap);

      n = vsnprintf(msg, SCIP_MAXSTRLEN, formatstr, ap);
      if( n < 0 )
         msg[SCIP_MAXSTRLEN-1] = '\0';
      else if( n >= SCIP_MAXSTRLEN )
      {
         char* bigmsg;
#ifndef NDEBUG
         int m;
#endif
         
         if( BMSallocMemorySize(&bigmsg, n+1) == NULL )
         {
            va_end(aq);
            return;
         }
         
#ifndef NDEBUG
         m = vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#else
         vsnprintf(bigmsg, (size_t) n+1, formatstr, aq);
#endif
         assert(m == n);
         va_end(aq);
         messagePrintInfo(file, bigmsg, n);
         BMSfreeMemory(&bigmsg);
         return;
      }   
      messagePrintInfo(file, msg, SCIP_MAXSTRLEN);
      va_end(aq);
   }
}
