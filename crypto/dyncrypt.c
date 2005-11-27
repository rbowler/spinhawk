/*----------------------------------------------------------------------------*/
/* file: dyncrypt.c                                                           */
/*                                                                            */
/* Implementation of the z/Architecture crypto instructions described in      */
/* SA22-7832-04: z/Architecture Principles of Operation within the Hercules   */
/* z/Architecture emulator. This file may only be used with and within the    */
/* Hercules emulator for non-commercial use!                                  */
/*                                                                            */
/*                              (c) Copyright Bernard van der Helm, 2003-2005 */
/*                              Noordwijkerhout, The Netherlands.             */
/*----------------------------------------------------------------------------*/

#include "hstdinc.h"

#ifndef _DYNCRYPT_C_
#define _DYNCRYPT_C_
#endif 

#ifndef _DYNCRYPT_DLL_
#define _DYNCRYPT_DLL_
#endif

#include "hercules.h"
#include "opcode.h"
#include "inline.h"
#include "aes.h"
#include "des.h"
#include "sha1.h"
#include "sha256.h"

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST

/*----------------------------------------------------------------------------*/
/* Debugging options                                                          */
/*----------------------------------------------------------------------------*/
//#define OPTION_KIMD_DEBUG
//#define OPTION_KLMD_DEBUG
//#define OPTION_KM_DEBUG
//#define OPTION_KMAC_DEBUG
//#define OPTION_KMC_DEBUG

/*----------------------------------------------------------------------------*/
/* General Purpose Register 0 macro's (GR0)                                   */
/*----------------------------------------------------------------------------*/
/* fc   : Function code                                                       */
/* m    : Modifier bit                                                        */
/*----------------------------------------------------------------------------*/
#define GR0_fc(regs)    ((regs)->GR_L(0) & 0x0000007F)
#define GR0_m(regs)     (((regs)->GR_L(0) & 0x00000080) ? TRUE : FALSE)

/*----------------------------------------------------------------------------*/
/* Bit strings for query functions                                            */
/*----------------------------------------------------------------------------*/
#define KIMD_BITS       { 0xe0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#define KLMD_BITS       { 0xe0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#define KM_BITS         { 0xf0, 0x00, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#define KMAC_BITS       { 0xf0, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };
#define KMC_BITS        { 0xf0, 0x00, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

/*----------------------------------------------------------------------------*/
/* Write bytes on one line                                                    */
/*----------------------------------------------------------------------------*/
#define LOGBYTE(s, v, x) \
{ \
  int i; \
  \
  logmsg("  " s " "); \
  for(i = 0; i < (x); i++) \
    logmsg("%02X", (v)[i]); \
  logmsg("\n"); \
}

/*----------------------------------------------------------------------------*/
/* Write bytes on multiple lines                                              */
/*----------------------------------------------------------------------------*/
#define LOGBYTE2(s, v, x, y) \
{ \
  int i; \
  int j; \
  \
  logmsg("  " s "\n"); \
  for(i = 0; i < (y); i++) \
  { \
    logmsg("      "); \
    for(j = 0; j < (x); j++) \
      logmsg("%02X", (v)[i * (x) + j]); \
    logmsg("\n"); \
  } \
}

/*----------------------------------------------------------------------------*/
/* CPU determined amount of data (processed in one go)                        */
/*----------------------------------------------------------------------------*/
#define PROCESS_MAX     4096

/*----------------------------------------------------------------------------*/
/* Used for printing debugging info                                           */
/*----------------------------------------------------------------------------*/
#define TRUEFALSE(boolean)  ((boolean) ? "True" : "False")

#ifndef __SHA1_COMPILE__
#define __SHA1_COMPILE__
/*----------------------------------------------------------------------------*/
/* Get the chaining vector for output processing                              */
/*----------------------------------------------------------------------------*/
static void sha1_getcv(sha1_context *ctx, BYTE icv[20])
{
  int i, j;

  for(i = 0, j = 0; i < 5; i++)
  { 
    icv[j++] = (ctx->state[i] & 0xff000000) >> 24;
    icv[j++] = (ctx->state[i] & 0x00ff0000) >> 16;
    icv[j++] = (ctx->state[i] & 0x0000ff00) >> 8;
    icv[j++] = (ctx->state[i] & 0x000000ff);
  }
}

/*----------------------------------------------------------------------------*/
/* Set the initial chaining value                                             */
/*----------------------------------------------------------------------------*/
static void sha1_seticv(sha1_context *ctx, BYTE icv[20])
{
  int i, j;

  for(i = 0, j = 0; i < 5; i++)
  {
    ctx->state[i] = icv[j++] << 24;
    ctx->state[i] |= icv[j++] << 16;
    ctx->state[i] |= icv[j++] << 8;
    ctx->state[i] |= icv[j++];
  }
}
#endif /* __SHA1_COMPILE__ */

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
#ifndef __SHA256_COMPILE__
#define __SHA256_COMPILE__
/*----------------------------------------------------------------------------*/
/* Get the chaining vector for output processing                              */
/*----------------------------------------------------------------------------*/
static void sha256_getcv(sha256_context *ctx, BYTE icv[32])
{
  int i, j;

  for(i = 0, j = 0; i < 8; i++)
  {
    icv[j++] = (ctx->state[i] & 0xff000000) >> 24;
    icv[j++] = (ctx->state[i] & 0x00ff0000) >> 16;
    icv[j++] = (ctx->state[i] & 0x0000ff00) >> 8;
    icv[j++] = (ctx->state[i] & 0x000000ff);
  }
}

/*----------------------------------------------------------------------------*/
/* Set the initial chaining value                                             */
/*----------------------------------------------------------------------------*/
static void sha256_seticv(sha256_context *ctx, BYTE icv[32])
{
  int i, j;

  for(i = 0, j = 0; i < 8; i++)
  {
    ctx->state[i] = icv[j++] << 24;
    ctx->state[i] |= icv[j++] << 16;
    ctx->state[i] |= icv[j++] << 8;
    ctx->state[i] |= icv[j++];
  }
}
#endif /* __SHA256_COMPILE__ */
#endif /* FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1 */

/*----------------------------------------------------------------------------*/
/* Needed functions from sha1.c and sha256.c.                                 */
/* We do our own counting and padding, we only need the hashing.              */
/*----------------------------------------------------------------------------*/
void sha1_process(sha1_context *ctx, BYTE data[64]);

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
void sha256_process(sha256_context *ctx, BYTE data[64]);
#endif

/*----------------------------------------------------------------------------*/
/* B93E Compute intermediate message digest (KIMD) FC 0                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kimd_query)(int r1, int r2, REGS *regs)
{
  BYTE kimd_bits[] = KIMD_BITS;

  UNREFERENCED(r1);
  UNREFERENCED(r2);

#ifdef OPTION_KIMD_DEBUG
  logmsg("  KIMD: function 0: query\n");
#endif

  /* Store the parameter block */
  ARCH_DEP(vstorec)(kimd_bits, 15, GR_A(1, regs), 1, regs);

  /* Set condition code 0 */
  regs->psw.cc = 0;
}

/*----------------------------------------------------------------------------*/
/* B93E Compute intermediate message digest (KIMD) FC 1                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kimd_sha_1)(int r1, int r2, REGS *regs)
{
  BYTE buffer[64];
  int crypted;
  BYTE cv[20];
  sha1_context context;

  UNREFERENCED(r1);

#ifdef OPTION_KIMD_DEBUG
  logmsg("  KIMD: function 1: sha-1\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 64 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 19, ACCTYPE_WRITE, regs);

  /* Fetch and set initial chaining value */
  ARCH_DEP(vfetchc)(cv, 19, GR_A(1, regs), 1, regs);
  sha1_seticv(&context, cv);

#ifdef OPTION_KIMD_DEBUG
  LOGBYTE("icv   :", cv, 20);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 64 < PROCESS_MAX)
  {
    /* Fetch and process a block of data */
    ARCH_DEP(vfetchc)(buffer, 63, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KIMD_DEBUG
    LOGBYTE2("  input :", buffer, 16, 4);
#endif

    sha1_process(&context, buffer);

    /* Store the output chaining value */
    sha1_getcv(&context, cv);
    ARCH_DEP(vstorec)(cv, 19, GR_A(1, regs), 1, regs);

#ifdef OPTION_KIMD_DEBUG
    LOGBYTE("ocv   :", cv, 20);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 64);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 64);

#ifdef OPTION_KIMD_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
/*----------------------------------------------------------------------------*/
/* B93E Compute intermediate message digest (KIMD) FC 2                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kimd_sha_256)(int r1, int r2, REGS *regs)
{
  BYTE buffer[64];
  int crypted;
  BYTE cv[32];
  sha256_context context;

  UNREFERENCED(r1);

#ifdef OPTION_KIMD_DEBUG
  logmsg("  KIMD: function 2: sha-256\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 64 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 31, ACCTYPE_WRITE, regs);

  /* Fetch and set initial chaining value */
  ARCH_DEP(vfetchc)(cv, 31, GR_A(1, regs), 1, regs);
  sha256_seticv(&context, cv);

#ifdef OPTION_KIMD_DEBUG
  LOGBYTE("icv   :", cv, 32);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 64 < PROCESS_MAX)
  {
    /* Fetch and process a block of data */
    ARCH_DEP(vfetchc)(buffer, 63, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KIMD_DEBUG
    LOGBYTE2("  input :", buffer, 16, 4);
#endif

    sha256_process(&context, buffer);

    /* Store the output chaining value */
    sha256_getcv(&context, cv);
    ARCH_DEP(vstorec)(cv, 31, GR_A(1, regs), 1, regs);

#ifdef OPTION_KIMD_DEBUG
    LOGBYTE("ocv   :", cv, 32);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 64);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 64);

#ifdef OPTION_KIMD_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}
#endif /* FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1 */

/*----------------------------------------------------------------------------*/
/* B93F Compute last message digest (KLMD) FC 0                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(klmd_query)(int r1, int r2, REGS *regs)
{
  BYTE klmd_bits[] = KLMD_BITS;

  UNREFERENCED(r1);
  UNREFERENCED(r2);

#ifdef OPTION_KLMD_DEBUG
  logmsg("  KLMD: function 0: query\n");
#endif

  /* Store the parameter block */
  ARCH_DEP(vstorec)(klmd_bits, 15, GR_A(1, regs), 1, regs);

  /* Set condition code 0 */
  regs->psw.cc = 0;
}

/*----------------------------------------------------------------------------*/
/* B93F Compute last message digest (KLMD) FC 1                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(klmd_sha_1)(int r1, int r2, REGS *regs)
{
  BYTE buffer[64];
  int crypted;
  BYTE cv[20];
  sha1_context context;
  int i;

  UNREFERENCED(r1);

#ifdef OPTION_KLMD_DEBUG
  logmsg("  KLMD: function 1: sha-1\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 19, ACCTYPE_WRITE, regs);

  /* Fetch and set initial chaining value */
  ARCH_DEP(vfetchc)(cv, 19, GR_A(1, regs), 1, regs);
  sha1_seticv(&context, cv);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("icv   :", cv, 20);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 64 < PROCESS_MAX)
  {
    /* Check for last block */
    if(GR_A(r2 + 1, regs) < 64)
      break;

    /* Fetch and process a block of data */
    ARCH_DEP(vfetchc)(buffer, 63, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE2("input :", buffer, 16, 4);
#endif

    sha1_process(&context, buffer);

    /* Store the output chaining value */
    sha1_getcv(&context, cv);
    ARCH_DEP(vstorec)(cv, 19, GR_A(1, regs), 1, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE("ocv   :", cv, 20);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 64);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 64);

#ifdef OPTION_KLMD_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

  }

  /* Check if cpu determined amount of data is processed */
  if(GR_A(r2 + 1, regs) >= 64)
  {
    regs->psw.cc = 3;
    return;
  }

  /* Fetch and process possible last block of data */
  if(GR_A(r2 + 1, regs))
  {
    ARCH_DEP(vfetchc)(buffer, GR_A(r2 + 1, regs) - 1, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE("input :", buffer, (int) GR_A(r2 + 1, regs));
#endif

  }

  /* Do the padding */
  i = GR_A(r2 + 1, regs);
  if(i > 55)
  {
    buffer[i++] = 0x80;
    while(i < 64)
      buffer[i++] = 0x00;
    sha1_process(&context, buffer);
    for(i = 0; i < 56; i++)
      buffer[i] = 0;
  }
  else
  {
    buffer[i++] = 0x80;
    while(i < 56)
      buffer[i++] = 0x00;
  }

  /* Fetch and set the message bit length */
  ARCH_DEP(vfetchc)(&buffer[56], 7, GR_A(1, regs) + 20, 1, regs);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("mbl   :", &buffer[56], 8);
#endif

  /* Calculate and store the message digest */
  sha1_process(&context, buffer);
  sha1_getcv(&context, cv);
  ARCH_DEP(vstorec)(cv, 19, GR_A(1, regs), 1, regs);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("md    :", cv, 20);
#endif

  /* Update registers */
  SET_GR_A(r2, regs, GR_A(r2,regs) + GR_A(r2 + 1, regs));
  SET_GR_A(r2 + 1, regs, 0);

#ifdef OPTION_KLMD_DEBUG
  logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
  logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

  /* Set condition code */
  regs->psw.cc = 0;
}

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
/*----------------------------------------------------------------------------*/
/* B93F Compute last message digest (KLMD) FC 2                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(klmd_sha_256)(int r1, int r2, REGS *regs)
{
  BYTE buffer[64];
  int crypted;
  BYTE cv[32];
  sha256_context context;
  int i;

  UNREFERENCED(r1);

#ifdef OPTION_KLMD_DEBUG
  logmsg("  KLMD: function 2: sha-256\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 31, ACCTYPE_WRITE, regs);

  /* Fetch and set initial chaining value */
  ARCH_DEP(vfetchc)(cv, 31, GR_A(1, regs), 1, regs);
  sha256_seticv(&context, cv);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("icv   :", cv, 32);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 64 < PROCESS_MAX)
  {
    /* Check for last block */
    if(GR_A(r2 + 1, regs) < 64)
      break;

    /* Fetch and process a block of data */
    ARCH_DEP(vfetchc)(buffer, 63, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE2("input :", buffer, 16, 4);
#endif

    sha256_process(&context, buffer);

    /* Store the output chaining value */
    sha256_getcv(&context, cv);
    ARCH_DEP(vstorec)(cv, 31, GR_A(1, regs), 1, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE("ocv   :", cv, 32);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 64);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 64);

#ifdef OPTION_KLMD_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

  }

  /* Check if cpu determined amount of data is processed */
  if(GR_A(r2 + 1, regs) >= 64)
  {
    regs->psw.cc = 3;
    return;
  }

  /* Fetch and process possible last block of data */
  if(GR_A(r2 + 1, regs))
  {
    ARCH_DEP(vfetchc)(buffer, GR_A(r2 + 1, regs) - 1, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KLMD_DEBUG
    LOGBYTE("input :", buffer, (int) GR_A(r2 + 1, regs));
#endif

  }

  /* Do the padding */
  i = GR_A(r2 + 1, regs);
  if(i > 55)
  {
    buffer[i++] = 0x80;
    while(i < 64)
      buffer[i++] = 0x00;
    sha256_process(&context, buffer);
    for(i = 0; i < 56; i++)
      buffer[i] = 0;
  }
  else
  {
    buffer[i++] = 0x80;
    while(i < 56)
      buffer[i++] = 0x00;
  }

  /* Fetch and set the message bit length */
  ARCH_DEP(vfetchc)(&buffer[56], 7, GR_A(1, regs) + 32, 1, regs);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("mbl   :", &buffer[56], 8);
#endif

  /* Calculate and store the message digest */
  sha256_process(&context, buffer);
  sha256_getcv(&context, cv);
  ARCH_DEP(vstorec)(cv, 31, GR_A(1, regs), 1, regs);

#ifdef OPTION_KLMD_DEBUG
  LOGBYTE("md    :", cv, 32);
#endif

  /* Update registers */
  SET_GR_A(r2, regs, GR_A(r2,regs) + GR_A(r2 + 1, regs));
  SET_GR_A(r2 + 1, regs, 0);

#ifdef OPTION_KLMD_DEBUG
  logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
  logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

  /* Set condition code */
  regs->psw.cc = 0;
}
#endif /* FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1 */

/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM) FC 0                                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(km_query)(int r1, int r2, REGS *regs)
{
  BYTE km_bits[] = KM_BITS;

  UNREFERENCED(r1);
  UNREFERENCED(r2);

#ifdef OPTION_KM_DEBUG
  logmsg("  KM: function 0: query\n");
#endif

  /* Store the parameter block */
  ARCH_DEP(vstorec)(km_bits, 15, GR_A(1, regs), 1, regs);

  /* Set condition code 0 and return */
  regs->psw.cc = 0;
}

/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM) FC 1                                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(km_dea)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  des_context context;
  BYTE k[8];

#ifdef OPTION_KM_DEBUG
  logmsg("  KM: function 1: dea\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Fetch and set the cryptographic key */
  ARCH_DEP(vfetchc)(k, 7, GR_A(1, regs), 1, regs);
  des_set_key(&context, k);

#ifdef OPTION_KM_DEBUG
  LOGBYTE("k     :", k, 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
      des_decrypt(&context, buffer, buffer); 
    else
      des_encrypt(&context, buffer, buffer);

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs, GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KM_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM) FC 2                                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(km_tdea_128)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  des3_context context;
  BYTE k[16];

#ifdef OPTION_KM_DEBUG
  logmsg("  KM: function 2: tdea-128\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 15, GR_A(1, regs), 1, regs);
  des3_set_2keys(&context, &k[0], &k[8]);

#ifdef OPTION_KM_DEBUG
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
      des3_decrypt(&context, buffer, buffer);
    else
      des3_encrypt(&context, buffer, buffer);

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs,GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KM_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM) FC 3                                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(km_tdea_192)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  des3_context context;
  BYTE k[24];

#ifdef OPTION_KM_DEBUG
  logmsg("  KM: function 3: tdea-192\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 23, GR_A(1, regs), 1, regs);
  des3_set_3keys(&context, &k[0], &k[8], &k[16]);

#ifdef OPTION_KM_DEBUG
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
  LOGBYTE("k3    :", &k[16], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
      des3_decrypt(&context, buffer, buffer);
    else
      des3_encrypt(&context, buffer, buffer);

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs,GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs, GR_A(r2+1, regs) - 8);

#ifdef OPTION_KM_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM) FC 18                                             */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(km_aes_128)(int r1, int r2, REGS *regs)
{
  BYTE buffer[16];
  int crypted;
  aes_context context;
  BYTE k[16];

#ifdef OPTION_KM_DEBUG
  logmsg("  KM: function 18: aes-128\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 16)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 15, GR_A(1, regs), 1, regs);
  aes_set_key(&context, k, 128);

#ifdef OPTION_KM_DEBUG
  LOGBYTE("k     :", &k[0], 16);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 16 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 15, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("input :", buffer, 16);
#endif

    /* Do the job */
    if(GR0_m(regs))
      aes_decrypt(&context, buffer, buffer);
    else
      aes_encrypt(&context, buffer, buffer);

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 15, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KM_DEBUG
    LOGBYTE("output:", buffer, 16);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 16);
    if(r1 != r2)
      SET_GR_A(r2, regs,GR_A(r2,regs) + 16);
    SET_GR_A(r2 + 1, regs, GR_A(r2+1, regs) - 16);

#ifdef OPTION_KM_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}
#endif /* FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1 */

/*----------------------------------------------------------------------------*/
/* B91E Compute message authentication code (KMAC) FC 0                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmac_query)(int r1, int r2, REGS *regs)
{
  BYTE kmac_bits[] = KMAC_BITS;

  UNREFERENCED(r1);
  UNREFERENCED(r2);

#ifdef OPTION_KMAC_DEBUG
  logmsg("  KMAC: function 0: query\n");
#endif

  /* Store the parameter block */
  ARCH_DEP(vstorec)(kmac_bits, 15, GR_A(1, regs), 1, regs);

  /* Set condition code 0 */
  regs->psw.cc = 0;
}

/*----------------------------------------------------------------------------*/
/* B91E Compute message authentication code (KMAC) FC 1                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmac_dea)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  BYTE cv[8];
  int i;
  des_context context;
  BYTE k[8];

  UNREFERENCED(r1);

#ifdef OPTION_KMAC_DEBUG
  logmsg("  KMAC: function 1: dea\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic key */
  ARCH_DEP(vfetchc)(k, 7, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context, k);

#ifdef OPTION_KMAC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k     :", k, 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* XOR the message with chaining value */
    for(i = 0; i < 8; i++)
      buffer[i] ^= cv[i];

    /* Calculate the output chaining value */
    des_encrypt(&context, buffer, cv);

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(cv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("ocv   :", cv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs, GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KMAC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B91E Compute message authentication code (KMAC) FC 2                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmac_tdea_128)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  BYTE cv[8];
  des_context context1;
  des_context context2;
  int i;
  BYTE k[16];

  UNREFERENCED(r1);

#ifdef OPTION_KMAC_DEBUG
  logmsg("  KMAC: function 2: tdea-128\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 15, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context1, &k[0]);
  des_set_key(&context2, &k[8]);

#ifdef OPTION_KMAC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* XOR the message with chaining value */
    for(i = 0; i < 8; i++)
      buffer[i] ^= cv[i];  

    /* Calculate the output chaining value */
    des_encrypt(&context1, buffer, cv);
    des_decrypt(&context2, cv, cv);
    des_encrypt(&context1, cv, cv);

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(cv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("ocv   :", cv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KMAC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B91E Compute message authentication code (KMAC) FC 3                       */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmac_tdea_192)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  BYTE cv[8];
  des_context context1;
  des_context context2;
  des_context context3;
  int i;
  BYTE k[24];

  UNREFERENCED(r1);

#ifdef OPTION_KMAC_DEBUG
  logmsg("  KMAC: function 3: tdea-192\n");
#endif

  /* Check special conditions */
  if(!r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8 || GR0_m(regs))
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 23, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context1, &k[0]);
  des_set_key(&context2, &k[8]);
  des_set_key(&context3, &k[16]);

#ifdef OPTION_KMAC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
  LOGBYTE("k3    :", &k[16], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* XOR the message with chaining value */
    for(i = 0; i < 8; i++)
      buffer[i] ^= cv[i];

    /* Calculate the output chaining value */
    des_encrypt(&context1, buffer, cv);
    des_decrypt(&context2, cv, cv);
    des_encrypt(&context3, cv, cv);

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(cv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMAC_DEBUG
    LOGBYTE("ocv   :", cv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r2, regs,GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KMAC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 0                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_query)(int r1, int r2, REGS *regs)
{
  BYTE kmc_bits[] = KMC_BITS;

  UNREFERENCED(r1);
  UNREFERENCED(r2);

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 0: query\n");
#endif

  /* Store the parameter block */
  ARCH_DEP(vstorec)(kmc_bits, 15, GR_A(1, regs), 1, regs);

  /* Set condition code 0 */
  regs->psw.cc = 0;
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 1                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_dea)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  BYTE cv[8];
  BYTE ocv[8];
  des_context context;
  int i;
  BYTE k[8];

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 1: dea\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic key */
  ARCH_DEP(vfetchc)(k, 7, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context, k);

#ifdef OPTION_KMC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k     :", k, 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
    {

      /* Save the ocv */
      memcpy(ocv, buffer, 8);

      /* Decrypt and XOR */
      des_decrypt(&context, buffer, buffer);
      for(i = 0; i < 8; i++)
        buffer[i] ^= cv[i];
    }
    else
    {
      /* XOR and encrypt */
      for(i = 0; i < 8; i++)
        buffer[i] ^= cv[i];
      des_encrypt(&context, buffer, buffer);

      /* Save the ocv */
      memcpy(ocv, buffer, 8);
    }

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(ocv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("ocv   :", ocv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2,regs) + 8);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 8);

#ifdef OPTION_KMC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }

    /* Set cv for next 8 bytes */
    memcpy(cv, ocv, 8);
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 2                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_tdea_128)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  int i;
  BYTE cv[8];
  BYTE ocv[8];
  des_context context1;
  des_context context2;
  BYTE k[16];

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 2: tdea-128\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 15, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context1, &k[0]);
  des_set_key(&context2, &k[8]);

#ifdef OPTION_KMC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
    {

      /* Save the ocv */
      memcpy(ocv, buffer, 8);

      /* Decrypt and XOR */
      des_decrypt(&context1, buffer, buffer);
      des_encrypt(&context2, buffer, buffer);
      des_decrypt(&context1, buffer, buffer);
      for(i = 0; i < 8; i++)
        buffer[i] ^= cv[i];
    }
    else
    {
      /* XOR and encrypt */
      for(i = 0 ; i < 8; i++)
        buffer[i] ^= cv[i];
      des_encrypt(&context1, buffer, buffer);
      des_decrypt(&context2, buffer, buffer);
      des_encrypt(&context1, buffer, buffer);

      /* Save the ocv */
      memcpy(ocv, buffer, 8);
    }

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(ocv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("ocv   :", ocv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2, regs) + 8);
    SET_GR_A(r2 + 1, regs, GR_A(r2+1, regs) - 8);

#ifdef OPTION_KMC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }

    /* Set cv for next 8 bytes */
    memcpy(cv, ocv, 8);
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 3                               */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_tdea_192)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  int i;
  BYTE cv[8];
  BYTE ocv[8];
  des_context context1;
  des_context context2;
  des_context context3;
  BYTE k[24];

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 3: tdea-192\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 23, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context1, &k[0]);
  des_set_key(&context2, &k[8]);
  des_set_key(&context3, &k[16]);

#ifdef OPTION_KMC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
  LOGBYTE("k3    :", &k[16], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    if(GR0_m(regs))
    {

      /* Save the ocv */
      memcpy(ocv, buffer, 8);

      /* Decrypt and XOR */
      des_decrypt(&context3, buffer, buffer);
      des_encrypt(&context2, buffer, buffer);
      des_decrypt(&context1, buffer, buffer);
      for(i = 0; i < 8; i++)
        buffer[i] ^= cv[i];
    }
    else
    {
      /* XOR and encrypt */
      for(i = 0; i < 8; i++)
        buffer[i] ^= cv[i];
      des_encrypt(&context1, buffer, buffer);
      des_decrypt(&context2, buffer, buffer);
      des_encrypt(&context3, buffer, buffer);

      /* Save the ocv */
      memcpy(ocv, buffer, 8);
    }

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(ocv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("ocv   :", ocv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2, regs) + 8);
    SET_GR_A(r2 + 1, regs, GR_A(r2+1, regs) - 8);

#ifdef OPTION_KMC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }

    /* Set cv for next 8 bytes */
    memcpy(cv, ocv, 8);
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 18                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_aes_128)(int r1, int r2, REGS *regs)
{
  BYTE buffer[16];
  int crypted;
  BYTE cv[16];
  BYTE ocv[16];
  aes_context context;
  int i;
  BYTE k[16];

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 18: aes-128\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 16)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 15, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 15, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic key */
  ARCH_DEP(vfetchc)(k, 15, GR_A(1, regs) + 16, 1, regs);
  aes_set_key(&context, k, 128);

#ifdef OPTION_KMC_DEBUG
  LOGBYTE("icv   :", cv, 16);
  LOGBYTE("k     :", k, 16);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 16 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 15, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("input :", buffer, 16);
#endif

    /* Do the job */
    if(GR0_m(regs))
    {

      /* Save the ovc */
      memcpy(ocv, buffer, 16);

      /* Decrypt and XOR */
      aes_decrypt(&context, buffer, buffer);
      for(i = 0; i < 16; i++)
        buffer[i] ^= cv[i];
    }
    else
    {
      /* XOR and encrypt */
      for(i = 0; i < 16; i++)
        buffer[i] ^= cv[i];
      aes_encrypt(&context, buffer, buffer);

      /* Save the ovc */
      memcpy(ocv, buffer, 16);
    }

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 15, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("output:", buffer, 15);
#endif

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(ocv, 15, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("ocv   :", ocv, 15);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 16);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2,regs) + 16);
    SET_GR_A(r2 + 1, regs,GR_A(r2+1,regs) - 16);

#ifdef OPTION_KMC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }

    /* Set cv for next 8 bytes */
    memcpy(cv, ocv, 8);
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC) FC 67                              */
/*----------------------------------------------------------------------------*/
static void ARCH_DEP(kmc_prng)(int r1, int r2, REGS *regs)
{
  BYTE buffer[8];
  int crypted;
  int i;
  BYTE cv[8];
  BYTE ocv[8];
  BYTE tcv[8];
  des_context context1;
  des_context context2;
  des_context context3;
  BYTE k[24];

#ifdef OPTION_KMC_DEBUG
  logmsg("  KMC: function 67: prng\n");
#endif

  /* Check special conditions */
  if(!r1 || r1 & 0x01 || !r2 || r2 & 0x01 || GR_A(r2 + 1, regs) % 8)
    ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);

  /* Return with cc 0 on zero length */
  if(!GR_A(r2 + 1, regs))
  {
    regs->psw.cc = 0;
    return;
  }

  /* Test writeability output chaining value */
  ARCH_DEP(validate_operand)(GR_A(1,regs), 1, 7, ACCTYPE_WRITE, regs);

  /* Fetch the initial chaining value */
  ARCH_DEP(vfetchc)(cv, 7, GR_A(1, regs), 1, regs);

  /* Fetch and set the cryptographic keys */
  ARCH_DEP(vfetchc)(k, 23, GR_A(1, regs) + 8, 1, regs);
  des_set_key(&context1, &k[0]);
  des_set_key(&context2, &k[8]);
  des_set_key(&context3, &k[16]);

#ifdef OPTION_KMC_DEBUG
  LOGBYTE("icv   :", cv, 8);
  LOGBYTE("k1    :", &k[0], 8);
  LOGBYTE("k2    :", &k[8], 8);
  LOGBYTE("k3    :", &k[16], 8);
#endif

  /* Try to process the CPU-determined amount of data */
  crypted = 0;
  while(crypted += 8 < PROCESS_MAX)
  {
    /* Fetch a block of data */
    ARCH_DEP(vfetchc)(buffer, 7, GR_A(r2, regs), r2, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("input :", buffer, 8);
#endif

    /* Do the job */
    des_encrypt(&context1, buffer, buffer);
    des_decrypt(&context2, buffer, buffer);
    des_encrypt(&context3, buffer, buffer);

    /* Save the temporary cv */
    memcpy(tcv, buffer, 8);

    /* XOR */
    for(i = 0; i < 8; i++)
      buffer[i] ^= cv[i];

    des_encrypt(&context1, buffer, buffer);
    des_decrypt(&context2, buffer, buffer);
    des_encrypt(&context3, buffer, buffer);

    /* Store the output */
    ARCH_DEP(vstorec)(buffer, 7, GR_A(r1, regs), r1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("output:", buffer, 8);
#endif

    /* XOR */
    for(i = 0; i < 8; i++)
      buffer[i] ^= tcv[i];

    des_encrypt(&context1, buffer, buffer);
    des_decrypt(&context2, buffer, buffer);
    des_encrypt(&context3, buffer, buffer);

    /* Save the ocv */
    memcpy(ocv, buffer, 8);

    /* Store the output chaining value */
    ARCH_DEP(vstorec)(ocv, 7, GR_A(1, regs), 1, regs);

#ifdef OPTION_KMC_DEBUG
    LOGBYTE("ocv   :", ocv, 8);
#endif

    /* Update the registers */
    SET_GR_A(r1, regs,GR_A(r1,regs) + 8);
    if(r1 != r2)
      SET_GR_A(r2, regs, GR_A(r2, regs) + 8);
    SET_GR_A(r2 + 1, regs, GR_A(r2+1, regs) - 8);

#ifdef OPTION_KMC_DEBUG
    logmsg("  GR%02d  : " F_GREG "\n", r1, (regs)->GR(r1));
    logmsg("  GR%02d  : " F_GREG "\n", r2, (regs)->GR(r2));
    logmsg("  GR%02d  : " F_GREG "\n", r2 + 1, (regs)->GR(r2 + 1));
#endif

    /* check for end of data */
    if(!GR_A(r2 + 1, regs))
    {
      regs->psw.cc = 0;
      return;
    }

    /* Set cv for next 8 bytes */
    memcpy(cv, ocv, 8);
  }

  /* CPU-determined amount of data processed */
  regs->psw.cc = 3;
}
#endif /* FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1 */

/*----------------------------------------------------------------------------*/
/* B93E Compute intermediate message digest (KIMD)                            */
/*----------------------------------------------------------------------------*/
DEF_INST(compute_intermediate_message_digest_d)
{
  int r1;
  int r2;

  RRE(inst, regs, r1, r2);

#ifdef OPTION_KIMD_DEBUG
  logmsg("KIMD: compute intermediate message digest\n");
  logmsg("  r1        : GR%02d\n", r1);
  logmsg("    address : " F_VADR "\n", regs->GR(r1));
  logmsg("  r2        : GR%02d\n", r2);
  logmsg("    address : " F_VADR "\n", regs->GR(r2));
  logmsg("    length  : " F_GREG "\n", regs->GR(r2 + 1));
  logmsg("  GR00      : " F_GREG "\n", regs->GR(0));
  logmsg("    bit 56  : %s\n", TRUEFALSE(GR0_m(regs)));
  logmsg("    fc      : %d\n", GR0_fc(regs));
  logmsg("  GR01      : " F_GREG "\n", regs->GR(1));
#endif

  switch(GR0_fc(regs))
  {
    case 0:
      ARCH_DEP(kimd_query)(r1, r2, regs);
      break;

    case 1:
      ARCH_DEP(kimd_sha_1)(r1, r2, regs);
      break;

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
    case 2:
      ARCH_DEP(kimd_sha_256)(r1, r2, regs);
      break;
#endif

    default:
      ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);
      break;
  }
}

/*----------------------------------------------------------------------------*/
/* B93F Compute last message digest (KLMD)                                    */
/*----------------------------------------------------------------------------*/
DEF_INST(compute_last_message_digest_d)
{
  int r1;
  int r2;

  RRE(inst, regs, r1, r2);

#ifdef OPTION_KLMD_DEBUG
  logmsg("KLMD: compute last message digest\n");
  logmsg("  r1        : GR%02d\n", r1);
  logmsg("    address : " F_VADR "\n", regs->GR(r1));
  logmsg("  r2        : GR%02d\n", r2);
  logmsg("    address : " F_VADR "\n", regs->GR(r2));
  logmsg("    length  : " F_GREG "\n", regs->GR(r2 + 1));
  logmsg("  GR00      : " F_GREG "\n", regs->GR(0));
  logmsg("    bit 56  : %s\n", TRUEFALSE(GR0_m(regs)));
  logmsg("    fc      : %d\n", GR0_fc(regs));
  logmsg("  GR01      : " F_GREG "\n", regs->GR(1));
#endif

  switch(GR0_fc(regs))
  {
    case 0:
      ARCH_DEP(klmd_query)(r1, r2, regs);
      break;

    case 1:
      ARCH_DEP(klmd_sha_1)(r1, r2, regs);
      break;

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
    case 2:
      ARCH_DEP(klmd_sha_256)(r1, r2, regs);
      break;
#endif

    default:
      ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);
      break;
  }
}

/*----------------------------------------------------------------------------*/
/* B92E Cipher message (KM)                                                   */
/*----------------------------------------------------------------------------*/
DEF_INST(cipher_message_d)
{
  int r1;
  int r2;

  RRE(inst, regs, r1, r2);

#ifdef OPTION_KM_DEBUG
  logmsg("KM: cipher message\n");
  logmsg("  r1        : GR%02d\n", r1);
  logmsg("    address : " F_VADR "\n", regs->GR(r1));
  logmsg("  r2        : GR%02d\n", r2);
  logmsg("    address : " F_VADR "\n", regs->GR(r2));
  logmsg("    length  : " F_GREG "\n", regs->GR(r2 + 1));
  logmsg("  GR00      : " F_GREG "\n", regs->GR(0));
  logmsg("    m       : %s\n", TRUEFALSE(GR0_m(regs)));
  logmsg("    fc      : %d\n", GR0_fc(regs));
  logmsg("  GR01      : " F_GREG "\n", regs->GR(1));
#endif

  switch(GR0_fc(regs))
  {
    case 0:
      ARCH_DEP(km_query)(r1, r2, regs);
      break;

    case 1:
      ARCH_DEP(km_dea)(r1, r2, regs);
      break;

    case 2:
      ARCH_DEP(km_tdea_128)(r1, r2, regs);
      break;

    case 3:
      ARCH_DEP(km_tdea_192)(r1, r2, regs);
      break;

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
    case 18:
      ARCH_DEP(km_aes_128)(r1, r2, regs);
      break;
#endif

    default:
      ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);
      break;
  }
}

/*----------------------------------------------------------------------------*/
/* B91E Compute message authentication code (KMAC)                            */
/*----------------------------------------------------------------------------*/
DEF_INST(compute_message_authentication_code_d)
{
  int r1;
  int r2;

  RRE(inst, regs, r1, r2);

#ifdef OPTION_KMAC_DEBUG
  logmsg("KMAC: compute message authentication code\n");
  logmsg("  r2        : GR%02d\n", r2);
  logmsg("    address : " F_VADR "\n", regs->GR(r2));
  logmsg("    length  : " F_GREG "\n", regs->GR(r2 + 1));
  logmsg("  GR00      : " F_GREG "\n", regs->GR(0));
  logmsg("    bit 56  : %s\n", TRUEFALSE(GR0_m(regs)));
  logmsg("    fc      : %d\n", GR0_fc(regs));
  logmsg("  GR01      : " F_GREG "\n", regs->GR(1));
#endif

  switch(GR0_fc(regs))
  {
    case 0:
      ARCH_DEP(kmac_query)(r1, r2, regs);
      break;

    case 1:
      ARCH_DEP(kmac_dea)(r1, r2, regs);
      break;

    case 2:
      ARCH_DEP(kmac_tdea_128)(r1, r2, regs);
      break;

    case 3:
      ARCH_DEP(kmac_tdea_192)(r1, r2, regs);
      break;

    default:
      ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);
      break;
  }
}

/*----------------------------------------------------------------------------*/
/* B92F Cipher message with chaining (KMC)                                    */
/*----------------------------------------------------------------------------*/
DEF_INST(cipher_message_with_chaining_d)
{
  int r1;
  int r2;

  RRE(inst, regs, r1, r2);

#ifdef OPTION_KMC_DEBUG
  logmsg("KMC: cipher message with chaining\n");
  logmsg("  r1        : GR%02d\n", r1);
  logmsg("    address : " F_VADR "\n", regs->GR(r1));
  logmsg("  r2        : GR%02d\n", r2);
  logmsg("    address : " F_VADR "\n", regs->GR(r2));
  logmsg("    length  : " F_GREG "\n", regs->GR(r2 + 1));
  logmsg("  GR00      : " F_GREG "\n", regs->GR(0));
  logmsg("    m       : %s\n", TRUEFALSE(GR0_m(regs)));
  logmsg("    fc      : %d\n", GR0_fc(regs));
  logmsg("  GR01      : " F_GREG "\n", regs->GR(1));
#endif

  switch(GR0_fc(regs))
  {
    case 0:
      ARCH_DEP(kmc_query)(r1, r2, regs);
      break;

    case 1:
      ARCH_DEP(kmc_dea)(r1, r2, regs);
      break;

    case 2:
      ARCH_DEP(kmc_tdea_128)(r1, r2, regs);
      break;

    case 3:
      ARCH_DEP(kmc_tdea_192)(r1, r2, regs);
      break;

#ifdef FEATURE_MESSAGE_SECURITY_ASSIST_EXTENSION_1
    case 18:
      ARCH_DEP(kmc_aes_128)(r1, r2, regs);
      break;

    case 67:
      ARCH_DEP(kmc_prng)(r1, r2, regs);
      break;
#endif

    default:
      ARCH_DEP(program_interrupt)(regs, PGM_SPECIFICATION_EXCEPTION);
      break;
  }
}

#endif /*defined(FEATURE_MESSAGE_SECURITY_ASSIST)*/

#if !defined(_GEN_ARCH)

#if defined(_ARCHMODE2)
 #define  _GEN_ARCH _ARCHMODE2
 #include "dyncrypt.c"
#endif

#if defined(_ARCHMODE3)
 #undef   _GEN_ARCH
 #define  _GEN_ARCH _ARCHMODE3
 #include "dyncrypt.c"
#endif

HDL_DEPENDENCY_SECTION;
{
   HDL_DEPENDENCY (HERCULES);
   HDL_DEPENDENCY (REGS);
   HDL_DEPENDENCY (DEVBLK);
// HDL_DEPENDENCY (SYSBLK);
// HDL_DEPENDENCY (WEBBLK);
}
END_DEPENDENCY_SECTION;

HDL_REGISTER_SECTION;
{
#if defined(_390_FEATURE_MESSAGE_SECURITY_ASSIST)
  HDL_REGISTER(s390_cipher_message, s390_cipher_message_d);
  HDL_REGISTER(s390_cipher_message_with_chaining, s390_cipher_message_with_chaining_d);
  HDL_REGISTER(s390_compute_intermediate_message_digest, s390_compute_intermediate_message_digest_d);
  HDL_REGISTER(s390_compute_last_message_digest, s390_compute_last_message_digest_d);
  HDL_REGISTER(s390_compute_message_authentication_code, s390_compute_message_authentication_code_d);
#endif /*defined(_390_FEATURE_MESSAGE_SECURITY_ASSIST)*/

#if defined(_900_FEATURE_MESSAGE_SECURITY_ASSIST)
  HDL_REGISTER(z900_cipher_message, z900_cipher_message_d);
  HDL_REGISTER(z900_cipher_message_with_chaining, z900_cipher_message_with_chaining_d);
  HDL_REGISTER(z900_compute_intermediate_message_digest, z900_compute_intermediate_message_digest_d);
  HDL_REGISTER(z900_compute_last_message_digest, z900_compute_last_message_digest_d);
  HDL_REGISTER(z900_compute_message_authentication_code, z900_compute_message_authentication_code_d);
#endif /*defined(_900_FEATURE_MESSAGE_SECURITY_ASSIST)*/

  logmsg("Crypto module loaded (c) Copyright Bernard van der Helm, 2003-2005\n");

}
END_REGISTER_SECTION;

#endif /*!defined(_GEN_ARCH)*/
