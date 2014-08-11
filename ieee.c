/* IEEE.C       (c) Copyright Willem Konynenberg, 2001-2003          */
/*              (c) Copyright Roger Bowler and others, 2003-2013     */
/*              Hercules Binary (IEEE) Floating Point Instructions   */

/*-------------------------------------------------------------------*/
/* This module implements ESA/390 Binary Floating-Point (IEEE 754)   */
/* instructions as described in SA22-7201-05 ESA/390 Principles of   */
/* Operation and SA22-7832-08 z/Architecture Principles of Operation.*/
/*-------------------------------------------------------------------*/

/*
 * Hercules System/370, ESA/390, z/Architecture emulator
 * ieee.c
 * Binary (IEEE) Floating Point Instructions
 * Copyright (c) 2001-2009 Willem Konynenberg <wfk@xos.nl>
 * TCEB, TCDB and TCXB contributed by Per Jessen, 20 September 2001.
 * THDER,THDR by Roger Bowler, 19 July 2003.
 * Additional instructions by Roger Bowler, November 2004:
 *  LXDBR,LXDB,LXEBR,LXEB,LDXBR,LEXBR,CXFBR,CXGBR,CFXBR,CGXBR,
 *  MXDBR,MXDB,MDEBR,MDEB,MADBR,MADB,MAEBR,MAEB,MSDBR,MSDB,
 *  MSEBR,MSEB,DIEBR,DIDBR,TBEDR,TBDR.
 * Licensed under the Q Public License
 * For details, see html/herclic.html
 * Based very loosely on float.c by Peter Kuschnerus, (c) 2000-2006.
 * Converted to use J.R.Hauser's softfloat package - RB, Oct 2013.
 * Floating-point extension facility - RB, April 2014:
 *  CELFBR,CDLFBR,CXLFBR,CLFEBR,CLFDBR,CLFXBR,
 *  CELGBR,CDLGBR,CXLGBR,CLGEBR,CLGDBR,CLGXBR.
 */

#include "hstdinc.h"

#if !defined(_HENGINE_DLL_)
#define _HENGINE_DLL_
#endif

#if !defined(_IEEE_C_)
#define _IEEE_C_
#endif

#include "hercules.h"

#if defined(FEATURE_BINARY_FLOATING_POINT)

#include "opcode.h"
#include "inline.h"

/* Definitions of BFP rounding methods */
#define RM_DEFAULT_ROUNDING             0
#define RM_BIASED_ROUND_TO_NEAREST      1
#define RM_PREPARE_SHORTER_PRECISION    3
#define RM_ROUND_TO_NEAREST             4
#define RM_ROUND_TOWARD_ZERO            5
#define RM_ROUND_TOWARD_POS_INF         6
#define RM_ROUND_TOWARD_NEG_INF         7

/* Macro to generate program check if invalid BFP rounding method */
#undef BFPRM_CHECK
#if !defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
#define BFPRM_CHECK(x,regs) \
        {if (!((x)==0 || (x)==1 || ((x)>=4 && (x)<=7))) \
            {regs->program_interrupt(regs, PGM_SPECIFICATION_EXCEPTION);}}
#else /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/
#define BFPRM_CHECK(x,regs) \
        {if (!((x)==0 || (x)==1 || ((x)>=3 && (x)<=7))) \
            {regs->program_interrupt(regs, PGM_SPECIFICATION_EXCEPTION);}}
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

#if !defined(_IEEE_C)
/* Architecture independent code goes within this ifdef */
#include "milieu.h"
#include "softfloat.h"

/* Floating point exponent bias values */
#define FLOAT32_BIAS    127
#define FLOAT64_BIAS    1023
#define FLOAT128_BIAS   16383

#endif  /* !defined(_IEEE_C) */

/* externally defined architecture-dependent functions */
#define vfetch4 ARCH_DEP(vfetch4)
#define vfetch8 ARCH_DEP(vfetch8)

/* locally defined architecture-dependent functions */
#define float_exception_masked ARCH_DEP(float_exception_masked)
#define float_exception ARCH_DEP(float_exception)
#define vfetch_float32 ARCH_DEP(vfetch_float32)
#define vfetch_float64 ARCH_DEP(vfetch_float64)

#define add_ebfp ARCH_DEP(add_ebfp)
#define add_lbfp ARCH_DEP(add_lbfp)
#define add_sbfp ARCH_DEP(add_sbfp)
#define compare_ebfp ARCH_DEP(compare_ebfp)
#define compare_lbfp ARCH_DEP(compare_lbfp)
#define compare_sbfp ARCH_DEP(compare_sbfp)
#define divide_ebfp ARCH_DEP(divide_ebfp)
#define divide_lbfp ARCH_DEP(divide_lbfp)
#define divide_sbfp ARCH_DEP(divide_sbfp)
#define load_test_ebfp ARCH_DEP(load_test_ebfp)
#define load_test_lbfp ARCH_DEP(load_test_lbfp)
#define load_test_sbfp ARCH_DEP(load_test_sbfp)
#define load_neg_ebfp ARCH_DEP(load_neg_ebfp)
#define load_neg_lbfp ARCH_DEP(load_neg_lbfp)
#define load_neg_sbfp ARCH_DEP(load_neg_sbfp)
#define load_pos_ebfp ARCH_DEP(load_pos_ebfp)
#define load_pos_lbfp ARCH_DEP(load_pos_lbfp)
#define load_pos_sbfp ARCH_DEP(load_pos_sbfp)
#define multiply_ebfp ARCH_DEP(multiply_ebfp)
#define multiply_lbfp ARCH_DEP(multiply_lbfp)
#define multiply_sbfp ARCH_DEP(multiply_sbfp)
#define multiply_add_ebfp ARCH_DEP(multiply_add_ebfp)
#define multiply_add_lbfp ARCH_DEP(multiply_add_lbfp)
#define multiply_add_sbfp ARCH_DEP(multiply_add_sbfp)
#define multiply_subtract_ebfp ARCH_DEP(multiply_subtract_ebfp)
#define multiply_subtract_lbfp ARCH_DEP(multiply_subtract_lbfp)
#define multiply_subtract_sbfp ARCH_DEP(multiply_subtract_sbfp)
#define squareroot_ebfp ARCH_DEP(squareroot_ebfp)
#define squareroot_lbfp ARCH_DEP(squareroot_lbfp)
#define squareroot_sbfp ARCH_DEP(squareroot_sbfp)
#define subtract_ebfp ARCH_DEP(subtract_ebfp)
#define subtract_lbfp ARCH_DEP(subtract_lbfp)
#define subtract_sbfp ARCH_DEP(subtract_sbfp)
#define test_data_class_ebfp ARCH_DEP(test_data_class_ebfp)
#define test_data_class_lbfp ARCH_DEP(test_data_class_lbfp)
#define test_data_class_sbfp ARCH_DEP(test_data_class_sbfp)
#define divint_lbfp ARCH_DEP(divint_lbfp)
#define divint_sbfp ARCH_DEP(divint_sbfp)

/*
 * Convert from Softfloat IEEE exception to Pop IEEE exception
 * with suppression of exceptions according to mask bits:
 *  - mask bit 0x04 (XxC) suppresses the inexact exception
 */
static int float_exception_masked(REGS * regs, int mask)
{
    int fpc = 0;
    int dxc = 0;
    int exc = 0;
    int8 flags = float_get_exception_flags();

    if (flags & float_flag_inexact) {
        fpc = FPC_FLAG_SFX;
    }

    if (flags & float_flag_underflow) {
        fpc |= FPC_FLAG_SFU;
    } else if (flags & float_flag_overflow) {
        fpc |= FPC_FLAG_SFO;
    } else if (flags & float_flag_divbyzero) {
        fpc |= FPC_FLAG_SFZ;
    } else if (flags & float_flag_invalid) {
        fpc |= FPC_FLAG_SFI;
    }

    exc = (fpc & FPC_FLAG) & ((regs->fpc & FPC_MASK) >> 8);

#if defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
    /* Suppress inexact exception if the XxC mask bit is set */
    if (mask & 0x04) {
        exc &= ~FPC_FLAG_SFX;
    }
#else
    UNREFERENCED(mask);
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

    if (exc & FPC_FLAG_SFI) {
        dxc = DXC_IEEE_INVALID_OP;
    } else if (exc & FPC_FLAG_SFZ) {
        dxc = DXC_IEEE_DIV_ZERO;
    } else if (exc & FPC_FLAG_SFO) {
        dxc = (fpc & FPC_FLAG_SFX) ? DXC_IEEE_OF_INEX_TRUNC : DXC_IEEE_OF_EXACT;
    } else if (exc & FPC_FLAG_SFU) {
        dxc = (fpc & FPC_FLAG_SFX) ? DXC_IEEE_UF_INEX_TRUNC : DXC_IEEE_UF_EXACT;
    } else if (exc & FPC_FLAG_SFX) {
        dxc = DXC_IEEE_INEXACT_TRUNC;
    }

    if (exc) {
        regs->dxc = dxc;
        if (dxc == DXC_IEEE_DIV_ZERO || dxc == DXC_IEEE_INVALID_OP) {
            /* suppress operation */
            regs->program_interrupt(regs, PGM_DATA_EXCEPTION);
        }
        fpc &= ~exc;
        regs->fpc |= fpc;
        /*
         * Other operations need to take appropriate action
         * to complete the operation.
         */
        return PGM_DATA_EXCEPTION;
    } else {
        /* Set flags in FPC */
        regs->fpc |= fpc;
        /* have caller take default action */
        return 0;
    }

} /* end function float_exception_masked */

/*
 * Convert from Softfloat IEEE exception to Pop IEEE exception
 * without suppression of exceptions by mask bits
 */
static inline int float_exception(REGS * regs)
{
    return float_exception_masked(regs, 0);
} /* end function float_exception */

#if !defined(_IEEE_C)
/*
 * Set softfloat rounding mode according to BFP rounding mode mask
 */
void set_rounding_mode(U32 fpcreg, int mask)
{
    int brm;
    int8 flrm;

    /* If mask is zero, obtain rounding mode from FPC register */
    if (mask == RM_DEFAULT_ROUNDING)
        brm = ((fpcreg & FPC_BRM) >> FPC_BRM_SHIFT) + 4;
    else
        brm = mask;

    /* Convert BFP rounding mode to nearest equivalent FE rounding mode */
    switch (brm) {
    case RM_ROUND_TO_NEAREST: /* Round to nearest ties to even */
        flrm = float_round_nearest_even;
        break;
    case RM_ROUND_TOWARD_ZERO: /* Round toward zero */
        flrm = float_round_to_zero;
        break;
    case RM_ROUND_TOWARD_POS_INF: /* Round toward +infinity */
        flrm = float_round_up;
        break;
    case RM_ROUND_TOWARD_NEG_INF: /* Round toward -infinity */
        flrm = float_round_down;
        break;
    default:
        flrm = float_round_nearest_even;
        break;
    } /* end switch(brm) */

    /* Set rounding mode for softfloat */
    float_set_rounding_mode(flrm);

} /* end function set_rounding_mode */
#endif  /* !defined(_IEEE_C) */

#if !defined(_IEEE_C)
/*
 * Get softfloat operands from registers
 */
static inline void get_float32(float32 *op, U32 *fpr) {
    *op = *fpr;
}
static inline void get_float64(float64 *op, U32 *fpr) {
    *op = ((U64)fpr[0] << 32) | fpr[1];
}
static inline void get_float128(float128 *op, U32 *fpr) {
    op->high = ((U64)fpr[0] << 32) | fpr[1];
    op->low = ((U64)fpr[FPREX] << 32) | fpr[FPREX+1];
}
#endif  /* !defined(_IEEE_C) */

/*
 * Fetch softfloat operands from memory
 */
static inline void vfetch_float32(float32 *op, VADR addr, int arn, REGS *regs) {
    *op = vfetch4(addr, arn, regs);
}
static inline void vfetch_float64(float64 *op, VADR addr, int arn, REGS *regs) {
    *op = vfetch8(addr, arn, regs);
}

#if !defined(_IEEE_C)
/*
 * Put softfloat operands into registers
 */
static inline void put_float32(float32 *op, U32 *fpr) {
    *fpr = *op;
}
static inline void put_float64(float64 *op, U32 *fpr) {
    fpr[0] = (U32)(*op >> 32);
    fpr[1] = (U32)(*op & 0xFFFFFFFF);
}
static inline void put_float128(float128 *op, U32 *fpr) {
    fpr[0] = (U32)(op->high >> 32);
    fpr[1] = (U32)(op->high & 0xFFFFFFFF);
    fpr[FPREX] = (U32)(op->low >> 32);
    fpr[FPREX+1] = (U32)(op->low & 0xFFFFFFFF);
}
#define _IEEE_C
#endif  /* !defined(_IEEE_C) */

/*
 * Chapter 9. Floating-Point Overview and Support Instructions
 */

#if defined(FEATURE_FPS_EXTENSIONS)
#if !defined(_CBH_FUNC)
/*
 * Convert binary floating point to hexadecimal long floating point
 * save result into long register and return condition code
 * Roger Bowler, 19 July 2003
 */
static int cnvt_bfp_to_hfp (float64 *op, U32 *fpr)
{
    int exp;
    U64 fract;
    U32 r0, r1;
    int cc;
    int sign;

    sign = float64_is_neg(*op) ? 1 : 0;

    if (float64_is_nan(*op)) {
        r0 = 0x7FFFFFFF;
        r1 = 0xFFFFFFFF;
        cc = 3;
    } else if (float64_is_inf(*op)) {
        r0 = sign ? 0xFFFFFFFF : 0x7FFFFFFF;
        r1 = 0xFFFFFFFF;
        cc = 3;
    } else if (float64_is_zero(*op)) {
        r0 = sign ? 0x80000000 : 0;
        r1 = 0;
        cc = 0;
    } else if (float64_is_subnormal(*op)) {
        r0 = sign ? 0x80000000 : 0;
        r1 = 0;
        cc = sign ? 1 : 2;
    } else {
        //logmsg("ieee: exp=%d (X\'%3.3x\')\tfract=%16.16"I64_FMT"x\n",
        //        float64_exp(*op), float64_exp(*op), float64_fract(*op));
        /* Insert an implied 1. in front of the 52 bit binary
           fraction and lengthen the result to 56 bits */
        fract = (U64)(float64_fract(*op) | 0x10000000000000ULL) << 3;

        /* The binary exponent is equal to the biased exponent - 1023
           adjusted by 1 to move the point before the 56 bit fraction */
        exp = float64_exp(*op) - 1023 + 1;

        //logmsg("ieee: adjusted exp=%d\tfract=%16.16"I64_FMT"x\n", exp, fract);
        /* Shift the fraction right one bit at a time until
           the binary exponent becomes a multiple of 4 */
        while (exp & 3)
        {
            exp++;
            fract >>= 1;
        }
        //logmsg("ieee:  shifted exp=%d\tfract=%16.16"I64_FMT"x\n", exp, fract);

        /* Convert the binary exponent into a hexadecimal exponent
           by dropping the last two bits (which are now zero) */
        exp >>= 2;

        /* If the hexadecimal exponent is less than -64 then return
           a signed zero result with a non-zero condition code */
        if (exp < -64)
        {
            r0 = sign ? 0x80000000 : 0;
            r1 = 0;
            cc = sign ? 1 : 2;
        }
        /* If the hexadecimal exponent exceeds +63 then return
           a signed maximum result with condition code 3 */
        else if (exp > 63)
        {
            r0 = sign ? 0xFFFFFFFF : 0x7FFFFFFF;
            r1 = 0xFFFFFFFF;
            cc = 3;
        }
        else
        {
            /* Convert the hexadecimal exponent to a characteristic
               by adding 64 */
            exp += 64;

            /* Pack the exponent and the fraction into the result */
            r0 = (sign ? 1<<31 : 0) | (exp << 24) | (fract >> 32);
            r1 = fract & 0xFFFFFFFF;
            cc = sign ? 1 : 2;
        }
    }
    /* Store high and low halves of result into fp register array
       and return condition code */
    fpr[0] = r0;
    fpr[1] = r1;
    return cc;

} /* end function cnvt_bfp_to_hfp */

/*
 * Convert hexadecimal long floating point register to
 * binary floating point and return condition code
 * Roger Bowler, 28 Nov 2004
 */
static int cnvt_hfp_to_bfp (U32 *fpr, int rounding,
        int bfp_fractbits, int bfp_emax, int bfp_ebias,
        int *result_sign, int *result_exp, U64 *result_fract)
{
    BYTE sign;
    short expo;
    U64 fract;
    int roundup = 0;
    int cc;
    U64 b;

    /* Break the source operand into sign, characteristic, fraction */
    sign = fpr[0] >> 31;
    expo = (fpr[0] >> 24) & 0x007F;
    fract = ((U64)(fpr[0] & 0x00FFFFFF) << 32) | fpr[1];

    /* Determine whether to round up or down */
    switch (rounding) {
    case RM_BIASED_ROUND_TO_NEAREST:
    case RM_ROUND_TO_NEAREST: roundup = 0; break;
    case RM_DEFAULT_ROUNDING:
    case RM_ROUND_TOWARD_ZERO: roundup = 0; break;
    case RM_ROUND_TOWARD_POS_INF: roundup = (sign ? 0 : 1); break;
    case RM_ROUND_TOWARD_NEG_INF: roundup = sign; break;
    } /* end switch(rounding) */

    /* Convert HFP zero to BFP zero and return cond code 0 */
    if (fract == 0) /* a = -0 or +0 */
    {
        *result_sign = sign;
        *result_exp = 0;
        *result_fract = 0;
        return 0;
    }

    /* Set the condition code */
    cc = sign ? 1 : 2;

    /* Convert the HFP characteristic to a true binary exponent */
    expo = (expo - 64) * 4;

    /* Convert true binary exponent to a biased exponent */
    expo += bfp_ebias;

    /* Shift the fraction left until leftmost 1 is in bit 8 */
    while ((fract & 0x0080000000000000ULL) == 0)
    {
        fract <<= 1;
        expo -= 1;
    }

    /* Convert 56-bit fraction to 55-bit with implied 1 */
    expo--;
    fract &= 0x007FFFFFFFFFFFFFULL;

    if (expo < -(bfp_fractbits-1)) /* |a| < Dmin */
    {
        if (expo == -(bfp_fractbits-1) - 1)
        {
            if (rounding == RM_BIASED_ROUND_TO_NEAREST
                || rounding == RM_ROUND_TO_NEAREST)
                roundup = 1;
        }
        if (roundup) { expo = 0; fract = 1; } /* Dmin */
        else { expo = 0; fract = 0; } /* Zero */
    }
    else if (expo < 1) /* Dmin <= |a| < Nmin */
    {
        /* Reinstate implied 1 in preparation for denormalization */
        fract |= 0x0080000000000000ULL;

        /* Denormalize to get exponent back in range */
        fract >>= (expo + (bfp_fractbits-1));
        expo = 0;
    }
    else if (expo > (bfp_emax+bfp_ebias)) /* |a| > Nmax */
    {
        cc = 3;
        if (roundup) { /* Inf */
            expo = (bfp_emax+bfp_ebias) + 1;
            fract = 0;
        } else { /* Nmax */
            expo = (bfp_emax+bfp_ebias);
            fract = 0x007FFFFFFFFFFFFFULL - (((U64)1<<(1+(55-bfp_fractbits)))-1);
        } /* Nmax */
    } /* end Nmax < |a| */

    /* Set the result sign and exponent */
    *result_sign = sign;
    *result_exp = expo;

    /* Apply rounding before truncating to final fraction length */
    b = ( (U64)1 ) << ( 55 - bfp_fractbits);
    if (roundup && (fract & b))
    {
        fract += b;
    }

    /* Convert 55-bit fraction to result fraction length */
    *result_fract = fract >> (55-bfp_fractbits);

    return cc;
} /* end function cnvt_hfp_to_bfp */

#define _CBH_FUNC
#endif /*!defined(_CBH_FUNC)*/

/*-------------------------------------------------------------------*/
/* B359 THDR  - CONVERT BFP TO HFP (long)                      [RRE] */
/* Roger Bowler, 19 July 2003                                        */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_long_to_float_long_reg)
{
    int r1, r2;
    float64 op2;

    RRE(inst, regs, r1, r2);
    //logmsg("THDR r1=%d r2=%d\n", r1, r2);
    HFPREG2_CHECK(r1, r2, regs);

    /* Load long BFP operand from R2 register */
    get_float64(&op2, regs->fpr + FPR2I(r2));

    /* Convert to HFP register and set condition code */
    regs->psw.cc =
        cnvt_bfp_to_hfp (&op2,
                         regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_bfp_long_to_float_long_reg) */

/*-------------------------------------------------------------------*/
/* B358 THDER - CONVERT BFP TO HFP (short to long)             [RRE] */
/* Roger Bowler, 19 July 2003                                        */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_short_to_float_long_reg)
{
    int r1, r2;
    float32 op2;
    float64 lbfp_op2;

    RRE(inst, regs, r1, r2);
    //logmsg("THDER r1=%d r2=%d\n", r1, r2);
    HFPREG2_CHECK(r1, r2, regs);

    /* Load short BFP operand from R2 register */
    get_float32(&op2, regs->fpr + FPR2I(r2));

    /* Lengthen short BFP operand to long BFP */
    lbfp_op2 = float32_to_float64(op2);

    /* Convert to HFP register and set condition code */
    regs->psw.cc =
        cnvt_bfp_to_hfp (&lbfp_op2,
                         regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_bfp_short_to_float_long_reg) */

/*-------------------------------------------------------------------*/
/* B351 TBDR  - CONVERT HFP TO BFP (long)                      [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_float_long_to_bfp_long_reg)
{
    int r1, r2, m3;
    float64 op1;
    int sign, exp;
    U64 fract;

    RRF_M(inst, regs, r1, r2, m3);
    //logmsg("TBDR r1=%d r2=%d\n", r1, r2);
    HFPREG2_CHECK(r1, r2, regs);
    BFPRM_CHECK(m3,regs);

    regs->psw.cc =
        cnvt_hfp_to_bfp (regs->fpr + FPR2I(r2), m3,
            /*fractbits*/52, /*emax*/1023, /*ebias*/1023,
            &sign, &exp, &fract);
    op1 = float64_build(sign, exp, fract);

    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_float_long_to_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B350 TBEDR - CONVERT HFP TO BFP (long to short)             [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_float_long_to_bfp_short_reg)
{
    int r1, r2, m3;
    float32 op1;
    int sign, exp;
    U64 fract;

    RRF_M(inst, regs, r1, r2, m3);
    //logmsg("TBEDR r1=%d r2=%d\n", r1, r2);
    HFPREG2_CHECK(r1, r2, regs);
    BFPRM_CHECK(m3,regs);

    regs->psw.cc =
        cnvt_hfp_to_bfp (regs->fpr + FPR2I(r2), m3,
            /*fractbits*/23, /*emax*/127, /*ebias*/127,
            &sign, &exp, &fract);
    op1 = float32_build(sign, exp, (U32)fract);

    put_float32(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_float_long_to_bfp_short_reg) */
#endif /*defined(FEATURE_FPS_EXTENSIONS)*/

/*
 * Chapter 19. Binary-Floating-Point Instructions
 * Most of these instructions were defined as an update to ESA/390.
 * z/Architecture has added instructions for 64-bit integers.
 */

/*-------------------------------------------------------------------*/
/* ADD (extended)                                                    */
/*-------------------------------------------------------------------*/
static int add_ebfp(float128 *op1, float128 *op2, REGS *regs)
{
    int code;
    float128 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float128_add(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float128_is_nan(result) ? 3 :
                   float128_is_zero(result) ? 0 :
                   float128_is_neg(result) ? 1 : 2;
    return code;

} /* end function add_ebfp */

/*-------------------------------------------------------------------*/
/* B34A AXBR  - ADD (extended BFP)                             [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(add_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("AXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = add_ebfp(&op1, &op2, regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(add_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* ADD (long)                                                        */
/*-------------------------------------------------------------------*/
static int add_lbfp(float64 *op1, float64 *op2, REGS *regs)
{
    int code;
    float64 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float64_add(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float64_is_nan(result) ? 3 :
                   float64_is_zero(result) ? 0 :
                   float64_is_neg(result) ? 1 : 2;
    return code;

} /* end function add_lbfp */

/*-------------------------------------------------------------------*/
/* B31A ADBR  - ADD (long BFP)                                 [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(add_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("ADBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = add_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(add_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1A ADB   - ADD (long BFP)                                 [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(add_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("ADB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = add_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(add_bfp_long) */

/*-------------------------------------------------------------------*/
/* ADD (short)                                                       */
/*-------------------------------------------------------------------*/
static int add_sbfp(float32 *op1, float32 *op2, REGS *regs)
{
    int code;
    float32 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float32_add(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float32_is_nan(result) ? 3 :
                   float32_is_zero(result) ? 0 :
                   float32_is_neg(result) ? 1 : 2;
    return code;

} /* end function add_sbfp */

/*-------------------------------------------------------------------*/
/* B30A AEBR  - ADD (short BFP)                                [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(add_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("AEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = add_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(add_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED0A AEB   - ADD (short BFP)                                [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(add_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("AEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = add_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }
} /* end DEF_INST(add_bfp_short) */

/*-------------------------------------------------------------------*/
/* COMPARE (extended)                                                */
/*-------------------------------------------------------------------*/
static int compare_ebfp(float128 *op1, float128 *op2, int sig, REGS *regs)
{
    int code = 0;

    float_clear_exception_flags();
    if (float128_is_signaling_nan(*op1) || float128_is_signaling_nan(*op2)
        || (sig && (float128_is_nan(*op1) || float128_is_nan(*op2)))) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        if (code) {
            return code;
        }
    }

    regs->psw.cc = (float128_is_nan(*op1) || float128_is_nan(*op2)) ? 3 :
        float128_eq(*op1, *op2) ? 0 :
        float128_lt_quiet(*op1, *op2) ? 1 : 2;

    return code;

} /* end function compare_ebfp */

/*-------------------------------------------------------------------*/
/* B349 CXBR  - COMPARE (extended BFP)                         [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("CXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_ebfp(&op1, &op2, 0, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* COMPARE (long)                                                    */
/*-------------------------------------------------------------------*/
static int compare_lbfp(float64 *op1, float64 *op2, int sig, REGS *regs)
{
    int code = 0;

    float_clear_exception_flags();
    if (float64_is_signaling_nan(*op1) || float64_is_signaling_nan(*op2)
        || (sig && (float64_is_nan(*op1) || float64_is_nan(*op2)))) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        if (code) {
            return code;
        }
    }

    regs->psw.cc = (float64_is_nan(*op1) || float64_is_nan(*op2)) ? 3 :
        float64_eq(*op1, *op2) ? 0 :
        float64_lt_quiet(*op1, *op2) ? 1 : 2;

    return code;

} /* end function compare_lbfp */

/*-------------------------------------------------------------------*/
/* B319 CDBR  - COMPARE (long BFP)                             [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("CDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_lbfp(&op1, &op2, 0, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED19 CDB   - COMPARE (long BFP)                             [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("CDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = compare_lbfp(&op1, &op2, 0, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_bfp_long) */

/*-------------------------------------------------------------------*/
/* COMPARE (short)                                                   */
/*-------------------------------------------------------------------*/
static int compare_sbfp(float32 *op1, float32 *op2, int sig, REGS *regs)
{
    int code = 0;

    float_clear_exception_flags();
    if (float32_is_signaling_nan(*op1) || float32_is_signaling_nan(*op2)
        || (sig && (float32_is_nan(*op1) || float32_is_nan(*op2)))) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        if (code) {
            return code;
        }
    }

    regs->psw.cc = (float32_is_nan(*op1) || float32_is_nan(*op2)) ? 3 :
        float32_eq(*op1, *op2) ? 0 :
        float32_lt_quiet(*op1, *op2) ? 1 : 2;

    return code;

} /* end function compare_sbfp */

/*-------------------------------------------------------------------*/
/* B309 CEBR  - COMPARE (short BFP)                            [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("CEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_sbfp(&op1, &op2, 0, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED09 CEB   - COMPARE (short BFP)                            [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("CEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = compare_sbfp(&op1, &op2, 0, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_bfp_short) */

/*-------------------------------------------------------------------*/
/* B348 KXBR  - COMPARE AND SIGNAL (extended BFP)              [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_and_signal_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("KXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_ebfp(&op1, &op2, 1, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_and_signal_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B318 KDBR  - COMPARE AND SIGNAL (long BFP)                  [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_and_signal_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("KDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_lbfp(&op1, &op2, 1, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_and_signal_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED18 KDB   - COMPARE AND SIGNAL (long BFP)                  [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_and_signal_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("KDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = compare_lbfp(&op1, &op2, 1, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_and_signal_bfp_long) */

/*-------------------------------------------------------------------*/
/* B308 KEBR  - COMPARE AND SIGNAL (short BFP)                 [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_and_signal_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("KEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = compare_sbfp(&op1, &op2, 1, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_and_signal_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED08 KEB   - COMPARE AND SIGNAL (short BFP)                 [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(compare_and_signal_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("KEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = compare_sbfp(&op1, &op2, 1, regs);

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(compare_and_signal_bfp_short) */

/*-------------------------------------------------------------------*/
/* B396 CXFBR - CONVERT FROM FIXED (32 to extended BFP)        [RRE] */
/* B396 CXFBRA - CONVERT FROM FIXED (32 to extended BFP)       [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix32_to_bfp_ext_reg)
{
    int r1, r2, m3, m4;
    float128 op1;
    S32 op2;

    RRE_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CXFBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);
    BFPRM_CHECK(m3, regs);

    op2 = regs->GR_L(r2);
    op1 = int32_to_float128(op2);
    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_fix32_to_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B395 CDFBR - CONVERT FROM FIXED (32 to long BFP)            [RRE] */
/* B395 CDFBRA - CONVERT FROM FIXED (32 to long BFP)           [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix32_to_bfp_long_reg)
{
    int r1, r2, m3, m4;
    float64 op1;
    S32 op2;

    RRE_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CDFBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3, regs);

    op2 = regs->GR_L(r2);
    op1 = int32_to_float64(op2);
    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_fix32_to_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B394 CEFBR - CONVERT FROM FIXED (32 to short BFP)           [RRE] */
/* B394 CEFBRA - CONVERT FROM FIXED (32 to short BFP)          [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix32_to_bfp_short_reg)
{
    int r1, r2, m3, m4;
    float32 op1;
    S32 op2;
    int pgm_check;

    RRE_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CEFBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3, regs);

    op2 = regs->GR_L(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = int32_to_float32(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_fix32_to_bfp_short_reg) */

#if defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
/*-------------------------------------------------------------------*/
/* B392 CXLFBR - CONVERT FROM LOGICAL (32 to extended BFP)     [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u32_to_bfp_ext_reg)
{
    int r1, r2, m3, m4;
    float128 op1;
    U32 op2;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CXLFBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_L(r2);
    op1 = uint32_to_float128(op2);
    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_u32_to_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B391 CDLFBR - CONVERT FROM LOGICAL (32 to long BFP)         [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u32_to_bfp_long_reg)
{
    int r1, r2, m3, m4;
    float64 op1;
    U32 op2;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CDLFBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_L(r2);
    op1 = uint32_to_float64(op2);
    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_u32_to_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B390 CELFBR - CONVERT FROM LOGICAL (32 to short BFP)        [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u32_to_bfp_short_reg)
{
    int r1, r2, m3, m4;
    float32 op1;
    U32 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CELFBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_L(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = uint32_to_float32(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_u32_to_bfp_short_reg) */
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A6 CXGBR - CONVERT FROM FIXED (64 to extended BFP)        [RRE] */
/* B3A6 CXGBRA - CONVERT FROM FIXED (64 to extended BFP)       [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix64_to_bfp_ext_reg)
{
    int r1, r2, m3, m4;
    float128 op1;
    S64 op2;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CXGBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);
    op1 = int64_to_float128(op2);
    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_fix64_to_bfp_ext_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A5 CDGBR - CONVERT FROM FIXED (64 to long BFP)            [RRE] */
/* B3A5 CDGBRA - CONVERT FROM FIXED (64 to long BFP)           [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix64_to_bfp_long_reg)
{
    int r1, r2, m3, m4;
    float64 op1;
    S64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CDGBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = int64_to_float64(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_fix64_to_bfp_long_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A4 CEGBR - CONVERT FROM FIXED (64 to short BFP)           [RRE] */
/* B3A4 CEGBRA - CONVERT FROM FIXED (64 to short BFP)          [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_fix64_to_bfp_short_reg)
{
    int r1, r2, m3, m4;
    float32 op1;
    S64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CEGBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = int64_to_float32(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_fix64_to_bfp_short_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A2 CXLGBR - CONVERT FROM LOGICAL (64 to extended BFP)     [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u64_to_bfp_ext_reg)
{
    int r1, r2, m3, m4;
    float128 op1;
    U64 op2;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CXLGBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);
    op1 = uint64_to_float128(op2);
    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(convert_u64_to_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B3A1 CDLGBR - CONVERT FROM LOGICAL (64 to long BFP)         [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u64_to_bfp_long_reg)
{
    int r1, r2, m3, m4;
    float64 op1;
    U64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CDLGBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = uint64_to_float64(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_u64_to_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B3A0 CELGBR - CONVERT FROM LOGICAL (64 to short BFP)        [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_u64_to_bfp_short_reg)
{
    int r1, r2, m3, m4;
    float32 op1;
    U64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CELGBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    op2 = regs->GR_G(r2);

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = uint64_to_float32(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_u64_to_bfp_short_reg) */
#endif /*defined(FEATURE_ESAME)*/
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

/*-------------------------------------------------------------------*/
/* B39A CFXBR - CONVERT TO FIXED (extended BFP to 32)          [RRF] */
/* B39A CFXBRA - CONVERT TO FIXED (extended BFP to 32)         [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_ext_to_fix32_reg)
{
    int r1, r2, m3, m4;
    S32 op1;
    float128 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CFXBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r2, regs);
    BFPRM_CHECK(m3,regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float128_to_int32(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float128_is_zero(op2) ? 0 :
        float128_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_ext_to_fix32_reg) */

/*-------------------------------------------------------------------*/
/* B399 CFDBR - CONVERT TO FIXED (long BFP to 32)              [RRF] */
/* B399 CFDBRA - CONVERT TO FIXED (long BFP to 32)             [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_long_to_fix32_reg)
{
    int r1, r2, m3, m4;
    S32 op1;
    float64 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CFDBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float64_to_int32(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float64_is_zero(op2) ? 0 :
        float64_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_long_to_fix32_reg) */

/*-------------------------------------------------------------------*/
/* B398 CFEBR - CONVERT TO FIXED (short BFP to 32)             [RRF] */
/* B398 CFEBRA - CONVERT TO FIXED (short BFP to 32)            [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_short_to_fix32_reg)
{
    int r1, r2, m3, m4;
    S32 op1;
    float32 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CFEBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float32_to_int32(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float32_is_zero(op2) ? 0 :
        float32_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_short_to_fix32_reg) */

#if defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
/*-------------------------------------------------------------------*/
/* B39E CLFXBR - CONVERT TO LOGICAL (extended BFP to 32)       [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_ext_to_u32_reg)
{
    int r1, r2, m3, m4;
    U32 op1;
    float128 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLFXBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r2, regs);
    BFPRM_CHECK(m3,regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float128_to_uint32(op2);
    pgm_check = float_exception_masked(regs, m4);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float128_is_zero(op2) ? 0 :
        float128_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_ext_to_u32_reg) */

/*-------------------------------------------------------------------*/
/* B39D CLFDBR - CONVERT TO LOGICAL (long BFP to 32)           [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_long_to_u32_reg)
{
    int r1, r2, m3, m4;
    U32 op1;
    float64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLFDBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float64_to_uint32(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float64_is_zero(op2) ? 0 :
        float64_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_long_to_u32_reg) */

/*-------------------------------------------------------------------*/
/* B39C CLFEBR - CONVERT TO LOGICAL (short BFP to 32)          [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_short_to_u32_reg)
{
    int r1, r2, m3, m4;
    U32 op1;
    float32 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLFEBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float32_to_uint32(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_L(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float32_is_zero(op2) ? 0 :
        float32_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_short_to_u32_reg) */
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3AA CGXBR - CONVERT TO FIXED (extended BFP to 64)          [RRF] */
/* B3AA CGXBRA - CONVERT TO FIXED (extended BFP to 64)         [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_ext_to_fix64_reg)
{
    int r1, r2, m3, m4;
    S64 op1;
    float128 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CGXBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r2, regs);
    BFPRM_CHECK(m3,regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float128_to_int64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float128_is_zero(op2) ? 0 :
        float128_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_ext_to_fix64_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A9 CGDBR - CONVERT TO FIXED (long BFP to 64)              [RRF] */
/* B3A9 CGDBRA - CONVERT TO FIXED (long BFP to 64)             [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_long_to_fix64_reg)
{
    int r1, r2, m3, m4;
    S64 op1;
    float64 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CGDBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float64_to_int64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float64_is_zero(op2) ? 0 :
        float64_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_long_to_fix64_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3A8 CGEBR - CONVERT TO FIXED (short BFP to 64)             [RRF] */
/* B3A8 CGEBRA - CONVERT TO FIXED (short BFP to 64)            [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_short_to_fix64_reg)
{
    int r1, r2, m3, m4;
    S64 op1;
    float32 op2;
    int pgm_check;

    RRF_MMA(inst, regs, r1, r2, m3, m4);
    //logmsg("CGEBR(A) r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float32_to_int64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float32_is_zero(op2) ? 0 :
        float32_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_short_to_fix64_reg) */
#endif /*defined(FEATURE_ESAME)*/

#if defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)
#if defined(FEATURE_ESAME)
/*-------------------------------------------------------------------*/
/* B3AE CLGXBR - CONVERT TO LOGICAL (extended BFP to 64)       [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_ext_to_u64_reg)
{
    int r1, r2, m3, m4;
    U64 op1;
    float128 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLGXBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r2, regs);
    BFPRM_CHECK(m3,regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float128_to_uint64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float128_is_zero(op2) ? 0 :
        float128_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_ext_to_u64_reg) */

/*-------------------------------------------------------------------*/
/* B3AD CLGDBR - CONVERT TO LOGICAL (long BFP to 64)           [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_long_to_u64_reg)
{
    int r1, r2, m3, m4;
    U64 op1;
    float64 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLGDBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float64_to_uint64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float64_is_zero(op2) ? 0 :
        float64_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_long_to_u64_reg) */

/*-------------------------------------------------------------------*/
/* B3AC CLGEBR - CONVERT TO LOGICAL (short BFP to 64)          [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(convert_bfp_short_to_u64_reg)
{
    int r1, r2, m3, m4;
    U64 op1;
    float32 op2;
    int pgm_check;

    RRF_MM(inst, regs, r1, r2, m3, m4);
    //logmsg("CLGEBR r1=%d r2=%d m3=%d m4=%d\n", r1, r2, m3, m4);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);

    op1 = float32_to_uint64(op2);
    pgm_check = float_exception_masked(regs, m4);

    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    regs->GR_G(r1) = op1;
    regs->psw.cc =
        (float_get_exception_flags() & float_flag_invalid) ? 3 :
        float32_is_zero(op2) ? 0 :
        float32_is_neg(op2) ? 1 : 2;

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(convert_bfp_short_to_u64_reg) */
#endif /*defined(FEATURE_ESAME)*/
#endif /*defined(FEATURE_FLOATING_POINT_EXTENSION_FACILITY)*/

/*-------------------------------------------------------------------*/
/* DIVIDE (extended)                                                 */
/*-------------------------------------------------------------------*/
static int divide_ebfp(float128 *op1, float128 *op2, REGS *regs)
{
    int code;
    float128 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float128_div(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function divide_ebfp */

/*-------------------------------------------------------------------*/
/* B34D DXBR  - DIVIDE (extended BFP)                          [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("DXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = divide_ebfp(&op1, &op2, regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* DIVIDE (long)                                                     */
/*-------------------------------------------------------------------*/
static int divide_lbfp(float64 *op1, float64 *op2, REGS *regs)
{
    int code;
    float64 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float64_div(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function divide_lbfp */

/*-------------------------------------------------------------------*/
/* B31D DDBR  - DIVIDE (long BFP)                              [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("DDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = divide_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1D DDB   - DIVIDE (long BFP)                              [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("DDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = divide_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_bfp_long) */

/*-------------------------------------------------------------------*/
/* DIVIDE (short)                                                    */
/*-------------------------------------------------------------------*/
static int divide_sbfp(float32 *op1, float32 *op2, REGS *regs)
{
    int code;
    float32 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float32_div(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function divide_sbfp */

/*-------------------------------------------------------------------*/
/* B30D DEBR  - DIVIDE (short BFP)                             [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("DEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = divide_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED0D DEB   - DIVIDE (short BFP)                             [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("DEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = divide_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_bfp_short) */

/*-------------------------------------------------------------------*/
/* B342 LTXBR - LOAD AND TEST (extended BFP)                   [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_and_test_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check = 0;

    RRE(inst, regs, r1, r2);
    //logmsg("LTXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    if (float128_is_signaling_nan(op2)) {
        float_raise(float_flag_invalid);
        pgm_check = float_exception(regs);
        op1 = float128_snan_to_qnan(op2);
    } else {
        op1 = op2;
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

    regs->psw.cc = float128_is_nan(op1) ? 3 :
                   float128_is_zero(op1) ? 0 :
                   float128_is_neg(op1) ? 1 : 2;

    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_and_test_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B312 LTDBR - LOAD AND TEST (long BFP)                       [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_and_test_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check = 0;

    RRE(inst, regs, r1, r2);
    //logmsg("LTDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    if (float64_is_signaling_nan(op2)) {
        float_raise(float_flag_invalid);
        pgm_check = float_exception(regs);
        op1 = float64_snan_to_qnan(op2);
    } else {
        op1 = op2;
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

    regs->psw.cc = float64_is_nan(op1) ? 3 :
                   float64_is_zero(op1) ? 0 :
                   float64_is_neg(op1) ? 1 : 2;

    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_and_test_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B302 LTEBR - LOAD AND TEST (short BFP)                      [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_and_test_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check = 0;

    RRE(inst, regs, r1, r2);
    //logmsg("LTEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    if (float32_is_signaling_nan(op2)) {
        float_raise(float_flag_invalid);
        pgm_check = float_exception(regs);
        op1 = float32_snan_to_qnan(op2);
    } else {
        op1 = op2;
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

    regs->psw.cc = float32_is_nan(op1) ? 3 :
                   float32_is_zero(op1) ? 0 :
                   float32_is_neg(op1) ? 1 : 2;

    put_float32(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_and_test_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* B357 FIEBR - LOAD FP INTEGER (short BFP)                    [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(load_fp_int_bfp_short_reg)
{
    int r1, r2, m3, pgm_check;
    float32 op1, op2;

    RRF_M(inst, regs, r1, r2, m3);
    //logmsg("FIEBR r1=%d, r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = float32_round_to_int(op2);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    pgm_check = float_exception(regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_fp_int_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* B35F FIDBR - LOAD FP INTEGER (long BFP)                     [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(load_fp_int_bfp_long_reg)
{
    int r1, r2, m3, pgm_check;
    float64 op1, op2;

    RRF_M(inst, regs, r1, r2, m3);
    //logmsg("FIDBR r1=%d, r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPRM_CHECK(m3,regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = float64_round_to_int(op2);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    pgm_check = float_exception(regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_fp_int_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B347 FIXBR - LOAD FP INTEGER (extended BFP)                 [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(load_fp_int_bfp_ext_reg)
{
    int r1, r2, m3, pgm_check;
    float128 op1, op2;

    RRF_M(inst, regs, r1, r2, m3);
    //logmsg("FIXBR r1=%d, r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);
    BFPRM_CHECK(m3,regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, m3);
    op1 = float128_round_to_int(op2);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);

    pgm_check = float_exception(regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_fp_int_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B29D LFPC  - LOAD FPC                                         [S] */
/* This instruction is in module esame.c                             */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* B304 LDEBR - LOAD LENGTHENED (short to long BFP)            [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_short_to_long_reg)
{
    int r1, r2;
    float64 op1;
    float32 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LDEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    float_clear_exception_flags();
    get_float32(&op2, regs->fpr + FPR2I(r2));

    op1 = float32_to_float64(op2);
    pgm_check = float_exception(regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_short_to_long_reg) */

/*-------------------------------------------------------------------*/
/* ED04 LDEB  - LOAD LENGTHENED (short to long BFP)            [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_short_to_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1;
    float32 op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("LDEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    float_clear_exception_flags();
    vfetch_float32(&op2, effective_addr2, b2, regs);

    op1 = float32_to_float64(op2);
    pgm_check = float_exception(regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_short_to_long) */

/*-------------------------------------------------------------------*/
/* B305 LXDBR - LOAD LENGTHENED (long to extended BFP)         [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_long_to_ext_reg)
{
    int r1, r2;
    float128 op1;
    float64 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LXDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    float_clear_exception_flags();
    get_float64(&op2, regs->fpr + FPR2I(r2));

    op1 = float64_to_float128(op2);
    pgm_check = float_exception(regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_long_to_ext_reg) */

/*-------------------------------------------------------------------*/
/* ED05 LXDB  - LOAD LENGTHENED (long to extended BFP)         [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_long_to_ext)
{
    int r1, b2;
    VADR effective_addr2;
    float128 op1;
    float64 op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("LXDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    float_clear_exception_flags();
    vfetch_float64(&op2, effective_addr2, b2, regs);

    op1 = float64_to_float128(op2);
    pgm_check = float_exception(regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_long_to_ext) */

/*-------------------------------------------------------------------*/
/* B306 LXEBR - LOAD LENGTHENED (short to extended BFP)        [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_short_to_ext_reg)
{
    int r1, r2;
    float128 op1;
    float32 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LXEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    float_clear_exception_flags();
    get_float32(&op2, regs->fpr + FPR2I(r2));

    op1 = float32_to_float128(op2);
    pgm_check = float_exception(regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_short_to_ext_reg) */

/*-------------------------------------------------------------------*/
/* ED06 LXEB  - LOAD LENGTHENED (short to extended BFP)        [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_lengthened_bfp_short_to_ext)
{
    int r1, b2;
    VADR effective_addr2;
    float128 op1;
    float32 op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("LXEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    float_clear_exception_flags();
    vfetch_float32(&op2, effective_addr2, b2, regs);

    op1 = float32_to_float128(op2);
    pgm_check = float_exception(regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_lengthened_bfp_short_to_ext) */

/*-------------------------------------------------------------------*/
/* B341 LNXBR - LOAD NEGATIVE (extended BFP)                   [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_negative_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LNXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    op1 = float128_neg(op2);

    regs->psw.cc = float128_is_nan(op1) ? 3 :
                   float128_is_zero(op1) ? 0 : 1;

    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_negative_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B311 LNDBR - LOAD NEGATIVE (long BFP)                       [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_negative_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LNDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    op1 = float64_neg(op2);

    regs->psw.cc = float64_is_nan(op1) ? 3 :
                   float64_is_zero(op1) ? 0 : 1;

    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_negative_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B301 LNEBR - LOAD NEGATIVE (short BFP)                      [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_negative_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LNEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    op1 = float32_neg(op2);

    regs->psw.cc = float32_is_nan(op1) ? 3 :
                   float32_is_zero(op1) ? 0 : 1;

    put_float32(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_negative_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* B343 LCXBR - LOAD COMPLEMENT (extended BFP)                 [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_complement_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LCXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    op1 = float128_is_neg(op2) ? float128_pos(op2) : float128_neg(op2);

    regs->psw.cc = float128_is_nan(op1) ? 3 :
                   float128_is_zero(op1) ? 0 :
                   float128_is_neg(op1) ? 1 : 2;

    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_complement_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B313 LCDBR - LOAD COMPLEMENT (long BFP)                     [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_complement_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LCDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    op1 = float64_is_neg(op2) ? float64_pos(op2) : float64_neg(op2);

    regs->psw.cc = float64_is_nan(op1) ? 3 :
                   float64_is_zero(op1) ? 0 :
                   float64_is_neg(op1) ? 1 : 2;

    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_complement_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B303 LCEBR - LOAD COMPLEMENT (short BFP)                    [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_complement_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LCEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    op1 = float32_is_neg(op2) ? float32_pos(op2) : float32_neg(op2);

    regs->psw.cc = float32_is_nan(op1) ? 3 :
                   float32_is_zero(op1) ? 0 :
                   float32_is_neg(op1) ? 1 : 2;

    put_float32(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_complement_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* B340 LPXBR - LOAD POSITIVE (extended BFP)                   [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_positive_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LPXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    op1 = float128_pos(op2);

    regs->psw.cc = float128_is_nan(op1) ? 3 :
                   float128_is_zero(op1) ? 0 : 2;

    put_float128(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_positive_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B310 LPDBR - LOAD POSITIVE (long BFP)                       [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_positive_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LPDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    op1 = float64_pos(op2);

    regs->psw.cc = float64_is_nan(op1) ? 3 :
                   float64_is_zero(op1) ? 0 : 2;

    put_float64(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_positive_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* B300 LPEBR - LOAD POSITIVE (short BFP)                      [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_positive_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;

    RRE(inst, regs, r1, r2);
    //logmsg("LPEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op2, regs->fpr + FPR2I(r2));

    op1 = float32_pos(op2);

    regs->psw.cc = float32_is_nan(op1) ? 3 :
                   float32_is_zero(op1) ? 0 : 2;

    put_float32(&op1, regs->fpr + FPR2I(r1));

} /* end DEF_INST(load_positive_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* B344 LEDBR - LOAD ROUNDED (long to short BFP)               [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_rounded_bfp_long_to_short_reg)
{
    int r1, r2;
    float32 op1;
    float64 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LEDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    op1 = float64_to_float32(op2);
    pgm_check = float_exception(regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check && (regs->fpc & (FPC_DXC_O | FPC_DXC_U))) {
        op2 = float32_to_float64(op1);
        put_float64(&op2, regs->fpr + FPR2I(r1));
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_rounded_bfp_long_to_short_reg) */

/*-------------------------------------------------------------------*/
/* B345 LDXBR - LOAD ROUNDED (extended to long BFP)            [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_rounded_bfp_ext_to_long_reg)
{
    int r1, r2;
    float64 op1;
    float128 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LDXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    op1 = float128_to_float64(op2);
    pgm_check = float_exception(regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check && (regs->fpc & (FPC_DXC_O | FPC_DXC_U))) {
        op2 = float64_to_float128(op1);
        put_float128(&op2, regs->fpr + FPR2I(r1));
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_rounded_bfp_ext_to_long_reg) */

/*-------------------------------------------------------------------*/
/* B346 LEXBR - LOAD ROUNDED (extended to short BFP)           [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(load_rounded_bfp_ext_to_short_reg)
{
    int r1, r2;
    float32 op1;
    float128 op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("LEXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op2, regs->fpr + FPR2I(r2));

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    op1 = float128_to_float32(op2);
    pgm_check = float_exception(regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check && (regs->fpc & (FPC_DXC_O | FPC_DXC_U))) {
        op2 = float32_to_float128(op1);
        put_float128(&op2, regs->fpr + FPR2I(r1));
    }

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(load_rounded_bfp_ext_to_short_reg) */

/*-------------------------------------------------------------------*/
/* MULTIPLY (extended)                                               */
/*-------------------------------------------------------------------*/
static int multiply_ebfp(float128 *op1, float128 *op2, REGS *regs)
{
    int code;
    float128 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float128_mul(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_ebfp */

/*-------------------------------------------------------------------*/
/* B34C MXBR  - MULTIPLY (extended BFP)                        [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("MXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = multiply_ebfp(&op1, &op2, regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* B307 MXDBR - MULTIPLY (long to extended BFP)                [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_long_to_ext_reg)
{
    int r1, r2;
    float64 op1, op2;
    float128 eb1, eb2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("MXDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    eb1 = float64_to_float128(op1);
    eb2 = float64_to_float128(op2);

    pgm_check = multiply_ebfp(&eb1, &eb2, regs);

    put_float128(&eb1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_long_to_ext_reg) */

/*-------------------------------------------------------------------*/
/* ED07 MXDB  - MULTIPLY (long to extended BFP)                [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_long_to_ext)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    float128 eb1, eb2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("MXDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    eb1 = float64_to_float128(op1);
    eb2 = float64_to_float128(op2);

    pgm_check = multiply_ebfp(&eb1, &eb2, regs);

    put_float128(&eb1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_long_to_ext) */

/*-------------------------------------------------------------------*/
/* MULTIPLY (long)                                                   */
/*-------------------------------------------------------------------*/
static int multiply_lbfp(float64 *op1, float64 *op2, REGS *regs)
{
    int code;
    float64 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float64_mul(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_lbfp */

/*-------------------------------------------------------------------*/
/* B31C MDBR  - MULTIPLY (long BFP)                            [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("MDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = multiply_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1C MDB   - MULTIPLY (long BFP)                            [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("MDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = multiply_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_long) */

/*-------------------------------------------------------------------*/
/* B30C MDEBR - MULTIPLY (short to long BFP)                   [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_short_to_long_reg)
{
    int r1, r2;
    float32 op1, op2;
    float64 lb1, lb2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("MDEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    lb1 = float32_to_float64(op1);
    lb2 = float32_to_float64(op2);

    pgm_check = multiply_lbfp(&lb1, &lb2, regs);

    put_float64(&lb1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_short_to_long_reg) */

/*-------------------------------------------------------------------*/
/* ED0C MDEB  - MULTIPLY (short to long BFP)                   [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_short_to_long)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    float64 lb1, lb2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("MDEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    lb1 = float32_to_float64(op1);
    lb2 = float32_to_float64(op2);

    pgm_check = multiply_lbfp(&lb1, &lb2, regs);

    put_float64(&lb1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_short_to_long) */

/*-------------------------------------------------------------------*/
/* MULTIPLY (short)                                                  */
/*-------------------------------------------------------------------*/
static int multiply_sbfp(float32 *op1, float32 *op2, REGS *regs)
{
    int code;
    float32 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float32_mul(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_sbfp */

/*-------------------------------------------------------------------*/
/* B317 MEEBR - MULTIPLY (short BFP)                           [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("MEEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = multiply_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED17 MEEB  - MULTIPLY (short BFP)                           [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("MEEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = multiply_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_bfp_short) */

/*-------------------------------------------------------------------*/
/* MULTIPLY AND ADD (long)                                           */
/*-------------------------------------------------------------------*/
static int multiply_add_lbfp(float64 *op1, float64 *op2, float64 *op3,
                REGS *regs)
{
    int code;
    float64 product, result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    product = float64_mul(*op2, *op3);
    result = float64_add(product, *op1);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_add_lbfp */

/*-------------------------------------------------------------------*/
/* B31E MADBR - MULTIPLY AND ADD (long BFP)                    [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_add_bfp_long_reg)
{
    int r1, r2, r3;
    float64 op1, op2, op3;
    int pgm_check;

    RRF_R(inst, regs, r1, r2, r3);
    //logmsg("MADBR r1=%d r3=%d r2=%d\n", r1, r3, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));
    get_float64(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_add_lbfp(&op1, &op2, &op3, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_add_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1E MADB  - MULTIPLY AND ADD (long BFP)                    [RXF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_add_bfp_long)
{
    int r1, r3, b2;
    VADR effective_addr2;
    float64 op1, op2, op3;
    int pgm_check;

    RXF(inst, regs, r1, r3, b2, effective_addr2);
    //logmsg("MADB r1=%d r3=%d b2=%d\n", r1, r3, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);
    get_float64(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_add_lbfp(&op1, &op2, &op3, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_add_bfp_long) */

/*-------------------------------------------------------------------*/
/* MULTIPLY AND ADD (short)                                          */
/*-------------------------------------------------------------------*/
static int multiply_add_sbfp(float32 *op1, float32 *op2, float32 *op3,
                REGS *regs)
{
    int code;
    float32 product, result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    product = float32_mul(*op2, *op3);
    result = float32_add(product, *op1);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_add_sbfp */

/*-------------------------------------------------------------------*/
/* B30E MAEBR - MULTIPLY AND ADD (short BFP)                   [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_add_bfp_short_reg)
{
    int r1, r2, r3;
    float32 op1, op2, op3;
    int pgm_check;

    RRF_R(inst, regs, r1, r2, r3);
    //logmsg("MAEBR r1=%d r3=%d r2=%d\n", r1, r3, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));
    get_float32(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_add_sbfp(&op1, &op2, &op3, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_add_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED0E MAEB  - MULTIPLY AND ADD (short BFP)                   [RXF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_add_bfp_short)
{
    int r1, r3, b2;
    VADR effective_addr2;
    float32 op1, op2, op3;
    int pgm_check;

    RXF(inst, regs, r1, r3, b2, effective_addr2);
    //logmsg("MAEB r1=%d r3=%d b2=%d\n", r1, r3, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);
    get_float32(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_add_sbfp(&op1, &op2, &op3, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_add_bfp_short) */

/*-------------------------------------------------------------------*/
/* MULTIPLY AND SUBTRACT (long)                                      */
/*-------------------------------------------------------------------*/
static int multiply_subtract_lbfp(float64 *op1, float64 *op2,
                float64 *op3, REGS *regs)
{
    int code;
    float64 product, result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    product = float64_mul(*op2, *op3);
    result = float64_sub(product, *op1);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_subtract_lbfp */

/*-------------------------------------------------------------------*/
/* B31F MSDBR - MULTIPLY AND SUBTRACT (long BFP)               [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_subtract_bfp_long_reg)
{
    int r1, r2, r3;
    float64 op1, op2, op3;
    int pgm_check;

    RRF_R(inst, regs, r1, r2, r3);
    //logmsg("MSDBR r1=%d r3=%d r2=%d\n", r1, r3, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));
    get_float64(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_subtract_lbfp(&op1, &op2, &op3, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_subtract_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1F MSDB  - MULTIPLY AND SUBTRACT (long BFP)               [RXF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_subtract_bfp_long)
{
    int r1, r3, b2;
    VADR effective_addr2;
    float64 op1, op2, op3;
    int pgm_check;

    RXF(inst, regs, r1, r3, b2, effective_addr2);
    //logmsg("MSDB r1=%d r3=%d b2=%d\n", r1, r3, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);
    get_float64(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_subtract_lbfp(&op1, &op2, &op3, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_subtract_bfp_long) */

/*-------------------------------------------------------------------*/
/* MULTIPLY AND SUBTRACT (short)                                     */
/*-------------------------------------------------------------------*/
static int multiply_subtract_sbfp(float32 *op1, float32 *op2,
                float32 *op3, REGS *regs)
{
    int code;
    float32 product, result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    product = float32_mul(*op2, *op3);
    result = float32_sub(product, *op1);
    code = float_exception(regs);
    *op1 = result;
    return code;

} /* end function multiply_subtract_sbfp */

/*-------------------------------------------------------------------*/
/* B30F MSEBR - MULTIPLY AND SUBTRACT (short BFP)              [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_subtract_bfp_short_reg)
{
    int r1, r2, r3;
    float32 op1, op2, op3;
    int pgm_check;

    RRF_R(inst, regs, r1, r2, r3);
    //logmsg("MSEBR r1=%d r3=%d r2=%d\n", r1, r3, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));
    get_float32(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_subtract_sbfp(&op1, &op2, &op3, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_subtract_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED0F MSEB  - MULTIPLY AND SUBTRACT (short BFP)              [RXF] */
/*-------------------------------------------------------------------*/
DEF_INST(multiply_subtract_bfp_short)
{
    int r1, r3, b2;
    VADR effective_addr2;
    float32 op1, op2, op3;
    int pgm_check;

    RXF(inst, regs, r1, r3, b2, effective_addr2);
    //logmsg("MSEB r1=%d r3=%d b2=%d\n", r1, r3, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);
    get_float32(&op3, regs->fpr + FPR2I(r3));

    pgm_check = multiply_subtract_sbfp(&op1, &op2, &op3, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(multiply_subtract_bfp_short) */

/*-------------------------------------------------------------------*/
/* B384 SFPC  - SET FPC                                        [RRE] */
/* This instruction is in module esame.c                             */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* B299 SRNM  - SET BFP ROUNDING MODE (2-BIT)                    [S] */
/* B2B8 SRNMB - SET BFP ROUNDING MODE (3-BIT)                    [S] */
/* These instructions are in module esame.c                          */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* SQUARE ROOT (extended)                                            */
/*-------------------------------------------------------------------*/
static int squareroot_ebfp(float128 *op, REGS *regs)
{
    int code;
    float128 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float128_sqrt(*op);
    code = float_exception(regs);
    *op = result;
    return code;

} /* end function squareroot_ebfp */

/*-------------------------------------------------------------------*/
/* B316 SQXBR - SQUARE ROOT (extended BFP)                     [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(squareroot_bfp_ext_reg)
{
    int r1, r2;
    float128 op;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SQXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op, regs->fpr + FPR2I(r2));

    pgm_check = squareroot_ebfp(&op, regs);

    put_float128(&op, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(squareroot_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* SQUARE ROOT (long)                                                */
/*-------------------------------------------------------------------*/
static int squareroot_lbfp(float64 *op, REGS *regs)
{
    int code;
    float64 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float64_sqrt(*op);
    code = float_exception(regs);
    *op = result;
    return code;

} /* end function squareroot_lbfp */

/*-------------------------------------------------------------------*/
/* B315 SQDBR - SQUARE ROOT (long BFP)                         [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(squareroot_bfp_long_reg)
{
    int r1, r2;
    float64 op;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SQDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op, regs->fpr + FPR2I(r2));

    pgm_check = squareroot_lbfp(&op, regs);

    put_float64(&op, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(squareroot_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED15 SQDB  - SQUARE ROOT (long BFP)                         [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(squareroot_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("SQDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    vfetch_float64(&op, effective_addr2, b2, regs);

    pgm_check = squareroot_lbfp(&op, regs);

    put_float64(&op, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(squareroot_bfp_long) */

/*-------------------------------------------------------------------*/
/* SQUARE ROOT (short)                                               */
/*-------------------------------------------------------------------*/
static int squareroot_sbfp(float32 *op, REGS *regs)
{
    int code;
    float32 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float32_sqrt(*op);
    code = float_exception(regs);
    *op = result;
    return code;

} /* end function squareroot_sbfp */

/*-------------------------------------------------------------------*/
/* B314 SQEBR - SQUARE ROOT (short BFP)                        [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(squareroot_bfp_short_reg)
{
    int r1, r2;
    float32 op;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SQEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op, regs->fpr + FPR2I(r2));

    pgm_check = squareroot_sbfp(&op, regs);

    put_float32(&op, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(squareroot_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED14 SQEB  - SQUARE ROOT (short BFP)                        [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(squareroot_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("SQEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    vfetch_float32(&op, effective_addr2, b2, regs);

    pgm_check = squareroot_sbfp(&op, regs);

    put_float32(&op, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(squareroot_bfp_short) */

/*-------------------------------------------------------------------*/
/* B29C STFPC - STORE FPC                                        [S] */
/* This instruction is in module esame.c                             */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* SUBTRACT (extended)                                               */
/*-------------------------------------------------------------------*/
static int subtract_ebfp(float128 *op1, float128 *op2, REGS *regs)
{
    int code;
    float128 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float128_sub(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float128_is_nan(result) ? 3 :
                   float128_is_zero(result) ? 0 :
                   float128_is_neg(result) ? 1 : 2;
    return code;

} /* end function subtract_ebfp */

/*-------------------------------------------------------------------*/
/* B34B SXBR  - SUBTRACT (extended BFP)                        [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(subtract_bfp_ext_reg)
{
    int r1, r2;
    float128 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SXBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);
    BFPREGPAIR2_CHECK(r1, r2, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));
    get_float128(&op2, regs->fpr + FPR2I(r2));

    pgm_check = subtract_ebfp(&op1, &op2, regs);

    put_float128(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(subtract_bfp_ext_reg) */

/*-------------------------------------------------------------------*/
/* SUBTRACT (long)                                                   */
/*-------------------------------------------------------------------*/
static int subtract_lbfp(float64 *op1, float64 *op2, REGS *regs)
{
    int code;
    float64 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float64_sub(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float64_is_nan(result) ? 3 :
                   float64_is_zero(result) ? 0 :
                   float64_is_neg(result) ? 1 : 2;
    return code;

} /* end function subtract_lbfp */

/*-------------------------------------------------------------------*/
/* B31B SDBR  - SUBTRACT (long BFP)                            [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(subtract_bfp_long_reg)
{
    int r1, r2;
    float64 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SDBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = subtract_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(subtract_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* ED1B SDB   - SUBTRACT (long BFP)                            [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(subtract_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("SDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    vfetch_float64(&op2, effective_addr2, b2, regs);

    pgm_check = subtract_lbfp(&op1, &op2, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(subtract_bfp_long) */

/*-------------------------------------------------------------------*/
/* SUBTRACT (short)                                                  */
/*-------------------------------------------------------------------*/
static int subtract_sbfp(float32 *op1, float32 *op2, REGS *regs)
{
    int code;
    float32 result;

    float_clear_exception_flags();
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    result = float32_sub(*op1, *op2);
    code = float_exception(regs);
    *op1 = result;
    regs->psw.cc = float32_is_nan(result) ? 3 :
                   float32_is_zero(result) ? 0 :
                   float32_is_neg(result) ? 1 : 2;
    return code;

} /* end function subtract_sbfp */

/*-------------------------------------------------------------------*/
/* B30B SEBR  - SUBTRACT (short BFP)                           [RRE] */
/*-------------------------------------------------------------------*/
DEF_INST(subtract_bfp_short_reg)
{
    int r1, r2;
    float32 op1, op2;
    int pgm_check;

    RRE(inst, regs, r1, r2);
    //logmsg("SEBR r1=%d r2=%d\n", r1, r2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = subtract_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(subtract_bfp_short_reg) */

/*-------------------------------------------------------------------*/
/* ED0B SEB   - SUBTRACT (short BFP)                           [RXE] */
/*-------------------------------------------------------------------*/
DEF_INST(subtract_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1, op2;
    int pgm_check;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("SEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    vfetch_float32(&op2, effective_addr2, b2, regs);

    pgm_check = subtract_sbfp(&op1, &op2, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(subtract_bfp_short) */

/*-------------------------------------------------------------------*/
/* ED10 TCEB  - TEST DATA CLASS (short BFP)                    [RXE] */
/* Per Jessen, Willem Konynenberg, 20 September 2001                 */
/*-------------------------------------------------------------------*/
DEF_INST(test_data_class_bfp_short)
{
    int r1, b2;
    VADR effective_addr2;
    float32 op1;
    int bit;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("TCEB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));

    if (float32_is_signaling_nan(op1))
        bit=30;
    else if (float32_is_nan(op1))
        bit=28;
    else if (float32_is_inf(op1))
        bit=26;
    else if (float32_is_subnormal(op1))
        bit=24;
    else if (float32_is_zero(op1))
        bit=20;
    else
        bit=22;

    if (float32_is_neg(op1))
        bit++;

    bit=31-bit;
    regs->psw.cc = (effective_addr2>>bit) & 1;

} /* end DEF_INST(test_data_class_bfp_short) */

/*-------------------------------------------------------------------*/
/* ED11 TCDB  - TEST DATA CLASS (long BFP)                     [RXE] */
/* Per Jessen, Willem Konynenberg, 20 September 2001                 */
/*-------------------------------------------------------------------*/
DEF_INST(test_data_class_bfp_long)
{
    int r1, b2;
    VADR effective_addr2;
    float64 op1;
    int bit;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("TCDB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));

    if (float64_is_signaling_nan(op1))
        bit=30;
    else if (float64_is_nan(op1))
        bit=28;
    else if (float64_is_inf(op1))
        bit=26;
    else if (float64_is_subnormal(op1))
        bit=24;
    else if (float64_is_zero(op1))
        bit=20;
    else
        bit=22;

    if (float64_is_neg(op1))
        bit++;

    bit=31-bit;
    regs->psw.cc = (effective_addr2>>bit) & 1;

} /* end DEF_INST(test_data_class_bfp_long) */

/*-------------------------------------------------------------------*/
/* ED12 TCXB  - TEST DATA CLASS (extended BFP)                 [RXE] */
/* Per Jessen, Willem Konynenberg, 20 September 2001                 */
/*-------------------------------------------------------------------*/
DEF_INST(test_data_class_bfp_ext)
{
    int r1, b2;
    VADR effective_addr2;
    float128 op1;
    int bit;

    RXE(inst, regs, r1, b2, effective_addr2);
    //logmsg("TCXB r1=%d b2=%d\n", r1, b2);
    BFPINST_CHECK(regs);
    BFPREGPAIR_CHECK(r1, regs);

    get_float128(&op1, regs->fpr + FPR2I(r1));

    if (float128_is_signaling_nan(op1))
        bit=30;
    else if (float128_is_nan(op1))
        bit=28;
    else if (float128_is_inf(op1))
        bit=26;
    else if (float128_is_subnormal(op1))
        bit=24;
    else if (float128_is_zero(op1))
        bit=20;
    else
        bit=22;

    if (float128_is_neg(op1))
        bit++;

    bit=31-bit;
    regs->psw.cc = (effective_addr2>>bit) & 1;

} /* end DEF_INST(test_data_class_bfp_ext) */

/*-------------------------------------------------------------------*/
/* DIVIDE TO INTEGER (long)                                          */
/*-------------------------------------------------------------------*/
static int divint_lbfp(float64 *op1, float64 *op2,
                        float64 *op3, int mode, REGS *regs)
{
    int code;
    int cc;
    int flags;
    int divsign, dvrsign;
    float128 xop1, xop2, xop3, xtemp;
    float128 xmaxint64 = {0x4034000000000000ULL, 0ULL}; // 2**53
    float128 xfract64 = {0xFFFFFFFFFFFFFFFFULL, 0xF000000000000000ULL}; // sign+exp+52 bits

    divsign = float64_is_neg(*op1) ? 1 : 0;
    dvrsign = float64_is_neg(*op2) ? 1 : 0;

    float_clear_exception_flags();

    if (float64_is_signaling_nan(*op1)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float64_snan_to_qnan(*op1);
        *op3 = float64_snan_to_qnan(*op1);
        regs->psw.cc = 1;
        return code;
    }

    if (float64_is_signaling_nan(*op2)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float64_snan_to_qnan(*op2);
        *op3 = float64_snan_to_qnan(*op2);
        regs->psw.cc = 1;
        return code;
    }

    if (float64_is_nan(*op1)) {
        *op3 = *op1;
        regs->psw.cc = 1;
        return 0;
    }

    if (float64_is_nan(*op2)) {
        *op1 = *op2;
        *op3 = *op2;
        regs->psw.cc = 1;
        return 0;
    }

    if (float64_is_inf(*op1) || float64_is_zero(*op2)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float64_default_nan;
        *op3 = float64_default_nan;
        regs->psw.cc = 1;
        return code;
    }

    if (float64_is_inf(*op2)) {
        *op3 = (divsign != dvrsign) ? 0x8000000000000000ULL : 0ULL;
        regs->psw.cc = 0;
        return 0;
    }

    xop1 = float64_to_float128(*op1);
    xop2 = float64_to_float128(*op2);

    set_rounding_mode(regs->fpc, RM_ROUND_TOWARD_ZERO);
    xtemp = float128_div(xop1, xop2);

    flags = float_get_exception_flags();
    if ((flags & float_flag_inexact)
        && (float128_le(xmaxint64, float128_pos(xtemp)))) {
        /* If quotient exceeds 2**53, truncate it to form a partial
           quotient consisting of an implied 1 with 52 fraction bits,
           and set the condition code to 2 */
        xop3.high = xtemp.high & xfract64.high;
        xop3.low = xtemp.low & xfract64.low;
        cc = 2;
    } else {
        /* Otherwise round it to an integer according to the specified
           rounding mode, and set the condition code to 0 */
        set_rounding_mode(regs->fpc, mode);
        float_clear_exception_flags();
        xop3 = float128_round_to_int(xtemp);
        cc = 0;
    }

    set_rounding_mode(regs->fpc, RM_ROUND_TOWARD_ZERO);
    float_clear_exception_flags();
    xtemp = float128_mul(xop2, xop3);
    float_clear_exception_flags();
    xop1 = float128_sub(xop1, xtemp);
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    float_clear_exception_flags();

    code = float_exception(regs);

    if (float128_exp(xop3) > FLOAT128_BIAS + FLOAT64_BIAS) {
        /* If quotient exponent exceeds maximum for float64, reduce
           the exponent by 1536 and increase condition code by 1 */
        xop3 = float128_build(float128_is_neg(xop3) ? 1 : 0,
                        float128_exp(xop3) - 1536,
                        float128_fract_high(xop3),
                        float128_fract_low(xop3));
        cc += 1;
    }

    *op1 = float128_to_float64(xop1);
    *op3 = float128_to_float64(xop3);

    /* A zero remainder is negative if the dividend is negative */
    if (float64_is_zero(*op1)) {
        *op1 = divsign ? 0x8000000000000000ULL : 0ULL;
    }

    /* A zero quotient is negative if dividend and divisor have different signs */
    if (float64_is_zero(*op3)) {
        *op3 = (divsign != dvrsign) ? 0x8000000000000000ULL : 0ULL;
    }

    regs->psw.cc = cc;
    return code;

} /* end function divint_lbfp */

/*-------------------------------------------------------------------*/
/* B35B DIDBR - DIVIDE TO INTEGER (long BFP)                   [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_integer_bfp_long_reg)
{
    int r1, r2, r3, m4;
    float64 op1, op2, op3;
    int pgm_check;

    RRF_RM(inst, regs, r1, r2, r3, m4);
    //logmsg("DIDBR r1=%d r3=%d r2=%d m4=%d\n", r1, r3, r2, m4);
    BFPINST_CHECK(regs);
    if (r1 == r2 || r2 == r3 || r1 == r3) {
        regs->program_interrupt(regs, PGM_SPECIFICATION_EXCEPTION);
    }
    BFPRM_CHECK(m4,regs);

    get_float64(&op1, regs->fpr + FPR2I(r1));
    get_float64(&op2, regs->fpr + FPR2I(r2));

    pgm_check = divint_lbfp(&op1, &op2, &op3, m4, regs);

    put_float64(&op1, regs->fpr + FPR2I(r1));
    put_float64(&op3, regs->fpr + FPR2I(r3));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_integer_bfp_long_reg) */

/*-------------------------------------------------------------------*/
/* DIVIDE TO INTEGER (short)                                         */
/*-------------------------------------------------------------------*/
static int divint_sbfp(float32 *op1, float32 *op2,
                        float32 *op3, int mode, REGS *regs)
{
    int code;
    int cc;
    int flags;
    int divsign, dvrsign;
    float128 xop1, xop2, xop3, xtemp;
    float128 xmaxint32 = {0x4017000000000000ULL, 0ULL}; // 2**24
    float128 xfract32 = {0xFFFFFFFFFE000000ULL, 0ULL}; // sign+exp+23 bits

    divsign = float32_is_neg(*op1) ? 1 : 0;
    dvrsign = float32_is_neg(*op2) ? 1 : 0;

    float_clear_exception_flags();

    if (float32_is_signaling_nan(*op1)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float32_snan_to_qnan(*op1);
        *op3 = float32_snan_to_qnan(*op1);
        regs->psw.cc = 1;
        return code;
    }

    if (float32_is_signaling_nan(*op2)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float32_snan_to_qnan(*op2);
        *op3 = float32_snan_to_qnan(*op2);
        regs->psw.cc = 1;
        return code;
    }

    if (float32_is_nan(*op1)) {
        *op3 = *op1;
        regs->psw.cc = 1;
        return 0;
    }

    if (float32_is_nan(*op2)) {
        *op1 = *op2;
        *op3 = *op2;
        regs->psw.cc = 1;
        return 0;
    }

    if (float32_is_inf(*op1) || float32_is_zero(*op2)) {
        float_raise(float_flag_invalid);
        code = float_exception(regs);
        *op1 = float32_default_nan;
        *op3 = float32_default_nan;
        regs->psw.cc = 1;
        return code;
    }

    if (float32_is_inf(*op2)) {
        *op3 = (divsign != dvrsign) ? 0x80000000 : 0;
        regs->psw.cc = 0;
        return 0;
    }

    xop1 = float32_to_float128(*op1);
    xop2 = float32_to_float128(*op2);

    set_rounding_mode(regs->fpc, RM_ROUND_TOWARD_ZERO);
    xtemp = float128_div(xop1, xop2);
    //logmsg("DIEBR div flags=%2.2X\n", float_get_exception_flags());

    flags = float_get_exception_flags();
    if ((flags & float_flag_inexact)
        && (float128_le(xmaxint32, float128_pos(xtemp)))) {
        /* If quotient exceeds 2**24, truncate it to form a partial
           quotient consisting of an implied 1 with 23 fraction bits,
           and set the condition code to 2 */
        xop3.high = xtemp.high & xfract32.high;
        xop3.low = 0;
        cc = 2;
    } else {
        /* Otherwise round it to an integer according to the specified
           rounding mode, and set the condition code to 0 */
        set_rounding_mode(regs->fpc, mode);
        float_clear_exception_flags();
        xop3 = float128_round_to_int(xtemp);
        //logmsg("DIEBR rou flags=%2.2X\n", float_get_exception_flags());
        cc = 0;
    }

    set_rounding_mode(regs->fpc, RM_ROUND_TOWARD_ZERO);
    float_clear_exception_flags();
    xtemp = float128_mul(xop2, xop3);
    //logmsg("DIEBR mul flags=%2.2X\n", float_get_exception_flags());
    float_clear_exception_flags();
    xop1 = float128_sub(xop1, xtemp);
    //logmsg("DIEBR sub flags=%2.2X\n", float_get_exception_flags());
    set_rounding_mode(regs->fpc, RM_DEFAULT_ROUNDING);
    float_clear_exception_flags();

    code = float_exception(regs);

    if (float128_exp(xop3) > FLOAT128_BIAS + FLOAT32_BIAS) {
        /* If quotient exponent exceeds maximum for float32, reduce
           the exponent by 192 and increase condition code by 1 */
        xop3 = float128_build(float128_is_neg(xop3) ? 1 : 0,
                        float128_exp(xop3) - 192,
                        float128_fract_high(xop3),
                        float128_fract_low(xop3));
        cc += 1;
    }

    *op1 = float128_to_float32(xop1);
    *op3 = float128_to_float32(xop3);

    /* A zero remainder is negative if the dividend is negative */
    if (float32_is_zero(*op1)) {
        *op1 = divsign ? 0x80000000 : 0;
    }

    /* A zero quotient is negative if dividend and divisor have different signs */
    if (float32_is_zero(*op3)) {
        *op3 = (divsign != dvrsign) ? 0x80000000 : 0;
    }

    regs->psw.cc = cc;
    return code;

} /* end function divint_sbfp */

/*-------------------------------------------------------------------*/
/* B353 DIEBR - DIVIDE TO INTEGER (short BFP)                  [RRF] */
/*-------------------------------------------------------------------*/
DEF_INST(divide_integer_bfp_short_reg)
{
    int r1, r2, r3, m4;
    float32 op1, op2, op3;
    int pgm_check;

    RRF_RM(inst, regs, r1, r2, r3, m4);
    //logmsg("DIEBR r1=%d r3=%d r2=%d m4=%d\n", r1, r3, r2, m4);
    BFPINST_CHECK(regs);
    if (r1 == r2 || r2 == r3 || r1 == r3) {
        regs->program_interrupt(regs, PGM_SPECIFICATION_EXCEPTION);
    }
    BFPRM_CHECK(m4,regs);

    get_float32(&op1, regs->fpr + FPR2I(r1));
    get_float32(&op2, regs->fpr + FPR2I(r2));

    pgm_check = divint_sbfp(&op1, &op2, &op3, m4, regs);

    put_float32(&op1, regs->fpr + FPR2I(r1));
    put_float32(&op3, regs->fpr + FPR2I(r3));

    if (pgm_check) {
        regs->program_interrupt(regs, pgm_check);
    }

} /* end DEF_INST(divide_integer_bfp_short_reg) */

#endif  /* FEATURE_BINARY_FLOATING_POINT */

#if !defined(_GEN_ARCH)

#if defined(_ARCHMODE2)
 #define  _GEN_ARCH _ARCHMODE2
 #include "ieee.c"
#endif

#if defined(_ARCHMODE3)
 #undef   _GEN_ARCH
 #define  _GEN_ARCH _ARCHMODE3
 #include "ieee.c"
#endif

#endif  /*!defined(_GEN_ARCH) */

/* end of ieee.c */
