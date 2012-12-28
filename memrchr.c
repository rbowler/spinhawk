/* MEMRCHR.C    (c) Copyright Volker Bandke, 2003                    */
/*              Hercules Right to Left memory scan                   */

/*-------------------------------------------------------------------*/
/*      Scans the memory block and reports the last occurrence of    */
/*      the specified byte in the buffer.  Returns a pointer to      */
/*      the byte if found, or NULL if not found.                     */
/*-------------------------------------------------------------------*/

#include "hstdinc.h"

#define _MEMRCHR_C_
#define _HUTIL_DLL_

#include "hercules.h"

#if !defined( HAVE_MEMRCHR )

#include "memrchr.h"

DLL_EXPORT void *memrchr(const void *buf, int c, size_t num)
{
   unsigned char *pMem;
   if (num == 0)
   {
      return NULL;
   }
   for (pMem = (unsigned char *) buf + num - 1; pMem >= (unsigned char *) buf; pMem--)
   {
      if (*pMem == (unsigned char) c) break;
   }
   if (pMem >= (unsigned char *) buf)
   {
      return ((void *) pMem);
   }
   return NULL;
}

#endif // !defined(HAVE_MEMRCHR)
