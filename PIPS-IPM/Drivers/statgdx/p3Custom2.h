/* $Id: p3Custom2.h 52165 2015-05-19 13:33:37Z sdirkse $
 * Header for C routines in p3Custom2.c
 * This code is sensitive to the order of the #defines, so inlining is not OK
 */

#ifndef _P3CUSTOM2_H
#define _P3CUSTOM2_H

/* all args (execName, libName, msg) are really short strings */
int xGetExecName (unsigned char *execName, unsigned char *msg);
int xGetLibName (unsigned char *libName, unsigned char *msg);

#endif /* ifndef _P3CUSTOM2_H */
