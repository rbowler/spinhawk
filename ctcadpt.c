/* CTCADPT.C    (c) Copyright James A. Pierson, 2002-2012            */
/*              (c) Copyright Roger Bowler, 2000-2012                */
/*              (c) Copyright Willem Konynenberg, 2000-2009          */
/*              (c) Copyright Vic Cross, 2001-2009                   */
/*              Hercules Channel-to-Channel Emulation Support        */

// vmnet     (C) Copyright Willem Konynenberg, 2000-2009
// CTCT      (C) Copyright Vic Cross, 2001-2009
// CTCT Enhanced (C) Copyright Peter J. Jansen, 2014

// Notes:
//   This module contains the remaining CTC emulation modes that
//   have not been moved to seperate modules. There is also logic
//   to allow old style 3088 device definitions for compatibility
//   and may be removed in a future release.
//
//   Please read README.NETWORKING for more info.
//

// Enhanced CTCT
// =============
//   Enhanced CTCT has been added here and not (yet) moved to a
//   separate module.  It is designed to support full 3088 CTCA
//   functionality (e.g. as needed for GRS and other components)
//   by communicating using a pair of TCP sockets with a likewise
//   configured Hercules instance on a different PC (or same PC).

//   The implementation is based mostly on an IBM publication,
//   "ESCON Channel-to-Channel Adapter", SA22-7203-00, although
//   no claim for completeness of this implemenation is feasible.

//   Enhanced CTCT functionality is configured using an unmodified
//   CTCT device statement syntax, but is selected when the MTU
//   bufsize >= CTCE_MTU_MIN (=32776) so that no new syntax is
//   introduced for this.  (The older CTCT type can still be used
//   as this typically used a much smaller MTU bufsize,
//   like 2048 or so.)
//   (Please note that 32776 = sizeof(CTCE_SOKPFX) + 32K, with
//   32K the maximum sCount experienced in CTC CCW programs.)

//   Enhanced CTCT requires an even-odd pair of port numbers per
//   device side but only the even port numbers are to be configured,
//   the odd port numbers are just derived by adding 1 to the
//   (configured) even port numbers.  The socket connection pairs
//   cross-connect, the arrows showing the send->receive direction :
//
//      x-lport-even -> y-rport-odd
//      x-lport-odd  <- y-rport-even
//
//   A sample CTCT device configuration is shown here:
//
//      Hercules PC Host A with IP address 192.168.1.100 :
//
//         0E40  3088  CTCT  30880  192.168.1.200  30880  32776
//         0E41  3088  CTCT  30882  192.168.1.200  30882  32776
//
//      Hercules PC Host B with IP address 192.168.1.200 :
//
//         0E40  3088  CTCT  30880  192.168.1.100  30880  32776
//         0E41  3088  CTCT  30882  192.168.1.100  30882  32776

#include "hstdinc.h"

#define _CTCADPT_C_
#define _HENGINE_DLL_

#include "hercules.h"
#include "devtype.h"
#include "ctcadpt.h"

#include "opcode.h"
#include "devtype.h"

// ====================================================================
// Constants
// ====================================================================

static char *CTCE_CmdStr[14] = {
    "PRE" , //  0 = 00 = Prepare
    "CTL" , //  1 = 01 = Control
    "RED" , //  2 = 02 = Read
    "WRT" , //  3 = 03 = Write
    "SCB" , //  4 = 04 = Sense Command Byte
    "???" , //  5 = 05 = Not Used
    "RBK" , //  6 = 06 = Read Backward
    "WEF" , //  7 = 07 = Write End Of File
    "NOP" , //  8 = 10 = No Operation
    "SEM" , //  9 = 11 = Set Extended Mode
    "SAS" , // 10 = 12 = Sense Adapter State
    "SID" , // 11 = 13 = Sense ID
    "RCD" , // 12 = 14 = Read Configuration Data
    "???"   // 13 = 15 = Invalid Command Code
};

static BYTE CTCE_Cmd[256] = {
    13, 3, 2, 8,10, 3, 2, 1,13, 3, 2, 8, 6, 3, 2, 1,
    13, 3, 2, 8, 4, 3, 2, 1,13, 3, 2, 8, 6, 3, 2, 1,
    13, 3, 2, 8,13, 3, 2, 1,13, 3, 2, 8, 6, 3, 2, 1,
    13, 3, 2, 8, 4, 3, 2, 1,13, 3, 2, 8, 6, 3, 2, 1,
    13, 3, 2,13,13, 3, 2, 1,13, 3, 2,13, 6, 3, 2, 1,
    13, 3, 2,13, 4, 3, 2, 1,13, 3, 2,13, 6, 3, 2, 1,
    13, 3, 2,13,13, 3, 2, 1,13, 3, 2,13, 6, 3, 2, 1,
    13, 3, 2,13, 4, 3, 2, 1,13, 3, 2,13, 6, 3, 2, 1,
    13, 7, 2, 8,13, 7, 2, 1,13, 7, 2, 8, 6, 7, 2, 1,
    13, 7, 2, 8, 4, 7, 2, 1,13, 7, 2, 8, 6, 7, 2, 1,
    13, 7, 2, 8,13, 7, 2, 1,13, 7, 2, 8, 6, 7, 2, 1,
    13, 7, 2, 8, 4, 7, 2, 1,13, 7, 2, 8, 6, 7, 2, 1,
    13, 7, 2, 9,13, 7, 2, 1,13, 7, 2,13, 6, 7, 2, 1,
    13, 7, 2,13, 4, 7, 2, 1,13, 7, 2,13, 6, 7, 2, 1,
    13, 7, 2, 0,11, 7, 2, 1,13, 7, 2,13, 6, 7, 2, 1,
    13, 7, 2,13, 4, 7, 2, 1,13, 7, 2,13, 6, 7, 2, 1
};

#define IS_CTCE_CCW_PRE(c)      ((CTCE_Cmd[c]==0))
#define IS_CTCE_CCW_CTL(c)      ((CTCE_Cmd[c]==1))
#define IS_CTCE_CCW_RED(c)      ((CTCE_Cmd[c]==2))
#define IS_CTCE_CCW_WRT(c)      ((CTCE_Cmd[c]==3))
#define IS_CTCE_CCW_SCB(c)      ((CTCE_Cmd[c]==4))
#define IS_CTCE_CCW_RBK(c)      ((CTCE_Cmd[c]==6))
#define IS_CTCE_CCW_WEF(c)      ((CTCE_Cmd[c]==7))
#define IS_CTCE_CCW_NOP(c)      ((CTCE_Cmd[c]==8))
#define IS_CTCE_CCW_SEM(c)      ((CTCE_Cmd[c]==9))
#define IS_CTCE_CCW_SAS(c)      ((CTCE_Cmd[c]==10))
#define IS_CTCE_CCW_SID(c)      ((CTCE_Cmd[c]==11))
#define IS_CTCE_CCW_RCD(c)      ((CTCE_Cmd[c]==12))
#define IS_CTCE_CCW_RDY(c)      ((CTCE_Cmd[c]<10))
#define IS_CTCE_CCW_RDA(c)      (((CTCE_Cmd[c]&0xFB)==2))

// ====================================================================
// Macros (might need relocation to a more appropriate .h file)
// ====================================================================

/* Macros for classifying CTC states          */
/* These are numbered 0 thru 7 as per the     */
/* column numbers 0-3 and 4-7 in the table    */
/* in section 2.13 in SA22-7203-00 by IBM     */
#define IS_CTCE_YWP(c)          (((c)&0x07)==0x00)
#define IS_CTCE_YWC(c)          (((c)&0x07)==0x01)
#define IS_CTCE_YWR(c)          (((c)&0x07)==0x02)
#define IS_CTCE_YWW(c)          (((c)&0x07)==0x03)
#define IS_CTCE_YAV(c)          (((c)&0x07)==0x04)
#define IS_CTCE_YNR(c)          (((c)&0x07)==0x05)
#define IS_CTCE_XWK(c)          (((c)&0x07)==0x06)
#define IS_CTCE_XIP(c)          (((c)&0x07)==0x07)

/* This last one is useful, tests for either  */
/* the 0 (YWP) or 4 (YAV)                     */
#define IS_CTCE_YAP(c)          (((c)&0x03)==0x00)

/* And the corresponding SET macros for these */
#define SET_CTCE_YWP(c)         (c&=&0x0F8)
#define SET_CTCE_YWC(c)         (c=(((c)&0xF8)|0x01))
#define SET_CTCE_YWR(c)         (c=(((c)&0xF8)|0x02))
#define SET_CTCE_YWW(c)         (c=(((c)&0xF8)|0x03))
#define SET_CTCE_YAV(c)         (c=(((c)&0xF8)|0x04))
#define SET_CTCE_YNR(c)         (c=(((c)&0xF8)|0x05))
#define SET_CTCE_XWK(c)         (c=(((c)&0xF8)|0x06))
#define SET_CTCE_XIP(c)         (c|=0x07))

/* Some additional flags are also present.    */
#define IS_CTCE_NRDY(c)         (((c)&0x80)==0x80)
#define IS_CTCE_WEOF(c)         (((c)&0x40)==0x40)
#define IS_CTCE_MATCH(c)        (((c)&0x20)==0x20)
#define IS_CTCE_ATTN(c)         (((c)&0x10)==0x10)

/* And the corresponding SET macros for these */
#define SET_CTCE_NRDY(c)        (c|=0x80)
#define SET_CTCE_WEOF(c)        (c|=0x40)
#define SET_CTCE_MATCH(c)       (c|=0x20)
#define SET_CTCE_ATTN(c)        (c|=0x10)

/* And the corresponding CLeaR macros         */
#define CLR_CTCE_NRDY(c)        (c&=~0x80)
#define CLR_CTCE_WEOF(c)        (c&=~0x40)
#define CLR_CTCE_MATCH(c)       (c&=~0x20)
#define CLR_CTCE_ATTN(c)        (c&=~0x10)

/* To CLeaR all flags                         */
#define CLR_CTCE_ALLF(c)        (c&=~0xF0)

/* The EOF and non-EOF Write CCW commands     */
#define IS_CCW_WRITE_EOF(c)     (((c)&0x83)==0x81)
#define IS_CCW_WRITE_NOEOF(c)   (((c)&0x83)==0x01)

/* Enhanced CTCT requires a minimum bufsize   */
/* in the Hercules configuration file, as the */
/* Send->Recv sockets may require this.  But  */
/* only when needed will the SeND packet be   */
/* larger than a SMaLl packet size.           */
#define CTCE_MTU_MIN            32776
#define CTCE_SND_SML            1500
//efine CTCE_SND_SML            32776

// ====================================================================
// Declarations
// ====================================================================

static int      CTCT_Init( DEVBLK *dev, int argc, char *argv[] );

static void     CTCT_Read( DEVBLK* pDEVBLK,   U16   sCount,
                           BYTE*   pIOBuf,    BYTE* pUnitStat,
                           U16*    pResidual, BYTE* pMore );

static void     CTCT_Write( DEVBLK* pDEVBLK,   U16   sCount,
                            BYTE*   pIOBuf,    BYTE* pUnitStat,
                            U16*    pResidual );

static void     CTCT_Send( DEVBLK* pDEVBLK,   U16   sCount,
                            BYTE*   pIOBuf,    BYTE* pUnitStat,
                            U16*    pResidual );

static void*    CTCT_RecvThread( void* argp );

static void*    CTCT_ListenThread( void* argp );

static int      VMNET_Init( DEVBLK *dev, int argc, char *argv[] );

static int      VMNET_Write( DEVBLK *dev, BYTE *iobuf,
                             U16 count, BYTE *unitstat );

static int      VMNET_Read( DEVBLK *dev, BYTE *iobuf,
                            U16 count, BYTE *unitstat );

// --------------------------------------------------------------------
// CTC Send-Receive Socket Prefix at the start of the DEVBLK buf
// --------------------------------------------------------------------

typedef struct _CTCE_SOKPFX
{
    BYTE                CmdReg;        /* CTCE command register      */
    BYTE                FsmSta;        /* CTCE FSM state             */
    U16                 sCount;        /* CTCE sCount copy           */
    U16                 PktCnt;        /* CTCE Packet Count          */
    U16                 SndLen;        /* CTCE Packet Sent Length    */
}
CTCE_SOKPFX;

// --------------------------------------------------------------------
// Definitions for CTC general data blocks
// --------------------------------------------------------------------

typedef struct _CTCG_PARMBLK
{
    int                 listenfd[2];   /* [0] = read, [1] = write    */
    u_int               ctceWrPort;    /*  0  = read,  1  = write    */
    struct sockaddr_in  addr;
    DEVBLK*             dev;
}
CTCG_PARMBLK;

// --------------------------------------------------------------------
// Device Handler Information Block
// --------------------------------------------------------------------

DEVHND ctcadpt_device_hndinfo =
{
        &CTCX_Init,                    /* Device Initialisation      */
        &CTCX_ExecuteCCW,              /* Device CCW execute         */
        &CTCX_Close,                   /* Device Close               */
        &CTCX_Query,                   /* Device Query               */
        NULL,                          /* Device Start channel pgm   */
        NULL,                          /* Device End channel pgm     */
        NULL,                          /* Device Resume channel pgm  */
        NULL,                          /* Device Suspend channel pgm */
        NULL,                          /* Device Read                */
        NULL,                          /* Device Write               */
        NULL,                          /* Device Query used          */
        NULL,                          /* Device Reserve             */
        NULL,                          /* Device Release             */
        NULL,                          /* Device Attention           */
        NULL,                          /* Immediate CCW Codes        */
        NULL,                          /* Signal Adapter Input       */
        NULL,                          /* Signal Adapter Output      */
        NULL,                          /* Hercules suspend           */
        NULL                           /* Hercules resume            */
};

DEVHND ctct_device_hndinfo =
{
        &CTCT_Init,                    /* Device Initialisation      */
        &CTCX_ExecuteCCW,              /* Device CCW execute         */
        &CTCX_Close,                   /* Device Close               */
        &CTCX_Query,                   /* Device Query               */
        NULL,                          /* Device Start channel pgm   */
        NULL,                          /* Device End channel pgm     */
        NULL,                          /* Device Resume channel pgm  */
        NULL,                          /* Device Suspend channel pgm */
        NULL,                          /* Device Read                */
        NULL,                          /* Device Write               */
        NULL,                          /* Device Query used          */
        NULL,                          /* Device Reserve             */
        NULL,                          /* Device Release             */
        NULL,                          /* Device Attention           */
        NULL,                          /* Immediate CCW Codes        */
        NULL,                          /* Signal Adapter Input       */
        NULL,                          /* Signal Adapter Output      */
        NULL,                          /* Hercules suspend           */
        NULL                           /* Hercules resume            */
};

DEVHND vmnet_device_hndinfo =
{
        &VMNET_Init,                   /* Device Initialisation      */
        &CTCX_ExecuteCCW,              /* Device CCW execute         */
        &CTCX_Close,                   /* Device Close               */
        &CTCX_Query,                   /* Device Query               */
        NULL,                          /* Device Start channel pgm   */
        NULL,                          /* Device End channel pgm     */
        NULL,                          /* Device Resume channel pgm  */
        NULL,                          /* Device Suspend channel pgm */
        NULL,                          /* Device Read                */
        NULL,                          /* Device Write               */
        NULL,                          /* Device Query used          */
        NULL,                          /* Device Reserve             */
        NULL,                          /* Device Release             */
        NULL,                          /* Device Attention           */
        NULL,                          /* Immediate CCW Codes        */
        NULL,                          /* Signal Adapter Input       */
        NULL,                          /* Signal Adapter Output      */
        NULL,                          /* Hercules suspend           */
        NULL                           /* Hercules resume            */
};

extern DEVHND ctci_device_hndinfo;
extern DEVHND lcs_device_hndinfo;

// ====================================================================
// Primary Module Entry Points
// ====================================================================

// --------------------------------------------------------------------
// Device Initialization Handler (Generic)
// --------------------------------------------------------------------

int  CTCX_Init( DEVBLK* pDEVBLK, int argc, char *argv[] )
{
    pDEVBLK->devtype = 0x3088;

    // The first argument is the device emulation type
    if( argc < 1 )
    {
        logmsg( _("HHCCT001E %4.4X: Incorrect number of parameters\n"),
            pDEVBLK->devnum );
        return -1;
    }

    if((pDEVBLK->hnd = hdl_ghnd(argv[0])))
    {
        if(pDEVBLK->hnd->init == &CTCX_Init)
            return -1;
        free(pDEVBLK->typname);
        pDEVBLK->typname = strdup(argv[0]);
        return (pDEVBLK->hnd->init)( pDEVBLK, --argc, ++argv );
    }
    logmsg (_("HHCCT034E %s: Unrecognized/unsupported CTC emulation type\n"),
        argv[0]);
    return -1;
}

// -------------------------------------------------------------------
// Query the device definition (Generic)
// -------------------------------------------------------------------

void  CTCX_Query( DEVBLK* pDEVBLK,
                  char**  ppszClass,
                  int     iBufLen,
                  char*   pBuffer )
{
    BEGIN_DEVICE_CLASS_QUERY( "CTCA", pDEVBLK, ppszClass, iBufLen, pBuffer );

    snprintf( pBuffer, iBufLen, "%s", pDEVBLK->filename );
}

// -------------------------------------------------------------------
// Close the device (Generic)
// -------------------------------------------------------------------

int  CTCX_Close( DEVBLK* pDEVBLK )
{
    // Close the device file (if not already closed)
    if( pDEVBLK->fd >= 0 )
    {
        if (socket_is_socket( pDEVBLK->fd ))
            close_socket( pDEVBLK->fd );
        else
            close( pDEVBLK->fd );
        pDEVBLK->fd = -1;           // indicate we're now closed
    }

    // Also for the RecvThread socket
    if( pDEVBLK->ctcemode && pDEVBLK->ctcefd >= 0 )
    {
        if (socket_is_socket( pDEVBLK->ctcefd ))
            close_socket( pDEVBLK->ctcefd );
        else
            close( pDEVBLK->ctcefd );
        pDEVBLK->ctcefd = -1;
    }

    return 0;
}

// -------------------------------------------------------------------
// Execute a Channel Command Word (Generic)
// -------------------------------------------------------------------

void  CTCX_ExecuteCCW( DEVBLK* pDEVBLK, BYTE  bCode,
                       BYTE    bFlags,  BYTE  bChained,
                       U16     sCount,  BYTE  bPrevCode,
                       int     iCCWSeq, BYTE* pIOBuf,
                       BYTE*   pMore,   BYTE* pUnitStat,
                       U16*    pResidual )
{
    int             iNum;               // Number of bytes to move
    BYTE            bOpCode;            // CCW opcode with modifier
                                        //   bits masked off

    UNREFERENCED( bFlags    );
    UNREFERENCED( bChained  );
    UNREFERENCED( bPrevCode );
    UNREFERENCED( iCCWSeq   );

    // Intervention required if the device file is not open
    if( ( pDEVBLK->fd < 0 || ( pDEVBLK->ctcemode && pDEVBLK->ctcefd < 0 ) ) &&
        !IS_CCW_SENSE( bCode ) &&
        !IS_CCW_CONTROL( bCode ) )
    {
        pDEVBLK->sense[0] = SENSE_IR;
        *pUnitStat = CSW_CE | CSW_DE | CSW_UC;
        return;
    }

    // Enhanced CTCT processing is implemented as an addition to
    // the already existing older code.
    if( pDEVBLK->ctcemode )
    {

        // Changes to DEVBLK are lock protected as the CTCT_RecvThread
        // might update as well, but, due to the way actually existing
        // software uses CTC devices, this may not be needed at all.
        obtain_lock( &pDEVBLK->lock );

        // Copy control command byte in x command register
        pDEVBLK->ctcexCmd = bCode;

        // The code below did not work, causing Wait State 064 RSN 9
        // at NIP time during IPL with GRS enabled in zOS.  Skipping
        // the initial NOT READY emulation (see CTCT_Init)  bypasses
        // the problem, and seemed not have any ill effects so far.
        // So this code remains disabled for the time being.
        //                                             (PJJ, April 2014)
#if 0
        // Commands different from sense adapter state (SAS) or sense ID
        // or read configuration data (RCD) with either or both sides of
        // the adapter NOT READY are processed first.
        if( ( IS_CTCE_NRDY( pDEVBLK->ctcexState ) ||
              IS_CTCE_NRDY( pDEVBLK->ctceyState ) ) &&
            IS_CTCE_CCW_RDY( bCode ) )
        {

            // If the x-side is not ready then we ensure it now will be
            if( IS_CTCE_NRDY( pDEVBLK->ctcexState ) )
            {
                CLR_CTCE_ALLF( pDEVBLK->ctcexState );
                SET_CTCE_YAV ( pDEVBLK->ctcexState );
            }

            // We always send the command to the other y-side where an
            // x-side not ready will generate a device end and reset
            // SENSE[0] bits 1 (SENSE_IR) and 7 (SENSE_OC) in CTCT_RecvThread.
            CTCT_Send( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );

            // If the y-side is not ready then unit check signals command reject.
            if( IS_CTCE_NRDY( pDEVBLK->ctceyState ) )
            {
                pDEVBLK->sense[0] = SENSE_CR;
                *pUnitStat = CSW_CE | CSW_DE | CSW_UC;
            }
            return;
        }
#endif

        // Any command different from READ and RDBACK will reset the WEOF flag.
//      if( !IS_CCW_READ( bCode ) &&
//          !IS_CCW_RDBACK( bCode ) )
        if( !IS_CTCE_CCW_RDA( bCode ) )
            CLR_CTCE_WEOF( pDEVBLK->ctcexState );

    }

    // Mask off the modifier bits in the CCW bOpCode
    if( ( bCode & 0x07 ) == 0x07 )
        bOpCode = 0x07;
    else if( ( bCode & 0x03 ) == 0x02 )
        bOpCode = 0x02;
    else if( ( bCode & 0x0F ) == 0x0C )
        bOpCode = 0x0C;
    else if( ( bCode & 0x03 ) == 0x01 )
        bOpCode = pDEVBLK->ctcxmode ? ( bCode & 0x83 ) : 0x01;
    else if( ( bCode & 0x1F ) == 0x14 )
        bOpCode = 0x14;
    else if( ( bCode & 0x47 ) == 0x03 )
        bOpCode = 0x03;
    else if( ( bCode & 0xC7 ) == 0x43 )
        bOpCode = 0x43;
    else
        bOpCode = bCode;

    // Process depending on CCW bOpCode
    switch (bOpCode)
    {
    case 0x01:  // 0MMMMM01  WRITE
        //------------------------------------------------------------
        // WRITE
        //------------------------------------------------------------

        // Enhanced CTCT processing.
        if( pDEVBLK->ctcemode )
        {

            // A write command is accepted if we are available and in
            // states YWR or YWP or YAV (the latter two equals YAP).
            if( IS_CTCE_YWR( pDEVBLK->ctcexState ) ||
                IS_CTCE_YAP( pDEVBLK->ctcexState ) )
            {

                // If available then we send the ATTN flag.
                if( IS_CTCE_YAV( pDEVBLK->ctcexState ) )
                {
                   SET_CTCE_ATTN( pDEVBLK->ctcexState );
                }

                // In state YWP or YAV (= YAP) then we move to state YWW.
                if( IS_CTCE_YAP( pDEVBLK->ctcexState ) )
                {
                    SET_CTCE_YWW( pDEVBLK->ctcexState );
                }

                // Otherwise we received a matching complementary command.
                else
                {
                    SET_CTCE_MATCH( pDEVBLK->ctcexState );
                }

                // Data and state info must be sent to the other (y-) side.
                CTCT_Send( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );

                // The above only proceeds when the matching read command
                // completes or was already done, unitistat and residual
                // will be set, so we just return to state available.
                SET_CTCE_YAV( pDEVBLK->ctcexState );
            }
            break;
        }
        // Non-Enhanced CTCT processing is retained below.

        // Return normal status if CCW count is zero
        if( sCount == 0 )
        {
            *pUnitStat = CSW_CE | CSW_DE;
            break;
        }

        // Write data and set unit status and residual byte count
        switch( pDEVBLK->ctctype )
        {
        case CTC_CTCT:
            CTCT_Write( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );
            break;
        case CTC_VMNET:
            *pResidual = sCount - VMNET_Write( pDEVBLK, pIOBuf,
                                               sCount,  pUnitStat );
            break;
        }
        break;

    case 0x81:  // 1MMMMM01  WEOF
        //------------------------------------------------------------
        // WRITE EOF
        //------------------------------------------------------------

        // Enhanced CTCT processing write EOF command is accepted if
        // in states YWR or YWP or YAV (the latter two equals YAP).
        if( ( pDEVBLK->ctcemode ) &&
            ( IS_CTCE_YWR( pDEVBLK->ctcexState ) ||
              IS_CTCE_YAP( pDEVBLK->ctcexState ) ) )
        {

            // This command is a matching complementary command for read.
            if( IS_CTCE_YWR( pDEVBLK->ctcexState ) )
                SET_CTCE_MATCH( pDEVBLK->ctcexState );

            // We then inform the other (y-)side we received a WEOF.
            CTCT_Send( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );

            // And we always return to state available.
            SET_CTCE_YAV( pDEVBLK->ctcexState );
            break;
        }
        // Non-Enhanced CTCT processing is retained.

        // Return normal status
        *pUnitStat = CSW_CE | CSW_DE;
        break;

    case 0x02:  // MMMMMM10  READ
    case 0x0C:  // MMMM1100  RDBACK
        // -----------------------------------------------------------
        // READ & READ BACKWARDS
        // -----------------------------------------------------------

        // Enhanced CTCT processing.
        if( pDEVBLK->ctcemode )
        {

            // If WEOF is set on our side whilst we are available and
            // the other side also available or in the W(D)P state, then
            // the read command will be rejected with unit exception.
            if( IS_CTCE_WEOF( pDEVBLK->ctcexState ) &&
                IS_CTCE_YAP ( pDEVBLK->ctcexState ) )
            {
                *pResidual = 0;
                *pUnitStat = CSW_CE | CSW_DE | CSW_UX ;
            }

            // A read command is accepted if we are available and in
            // states YWW or YWP or YAV (the latter two equals YAP).
            else if( IS_CTCE_YWW( pDEVBLK->ctcexState ) ||
                     IS_CTCE_YAP( pDEVBLK->ctcexState ) )
            {

                // If available then we send the ATTN flag.
                if( IS_CTCE_YAV( pDEVBLK->ctcexState ) )
                {
                    SET_CTCE_ATTN( pDEVBLK->ctcexState );
                }

                // In state YWP or YAV (= YAP) then we move to state YWR.
                if( IS_CTCE_YAP( pDEVBLK->ctcexState ) )
                {
                    SET_CTCE_YWR( pDEVBLK->ctcexState );
                }

                // Otherwise we received a matching complementary command.
                else
                {
                    SET_CTCE_MATCH( pDEVBLK->ctcexState );
                }

                // Data and state info must be sent to the other (y-) side,
                // proceeding only when the matching write command completes
                // or was already done, unitstat and residual will be set.
                CTCT_Send( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );
            }

            // And we return to state available.
            SET_CTCE_YAV( pDEVBLK->ctcexState );
            break;
        }
        // Non-Enhanced CTCT processing is retained below.

        // Read data and set unit status and residual byte count
        switch( pDEVBLK->ctctype )
        {
        case CTC_CTCT:
            CTCT_Read( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual, pMore );
            break;
        case CTC_VMNET:
            *pResidual = sCount - VMNET_Read( pDEVBLK, pIOBuf,
                                              sCount,  pUnitStat );
            break;
        }
        break;

    case 0x07:  // MMMMM111  CTL
        // -----------------------------------------------------------
        // CONTROL
        // -----------------------------------------------------------

        // Enhanced CTCT processing.
        if( pDEVBLK->ctcemode )
        {

            // The control command is accepted if we are available and
            // the other side is also available or in the YWP state.
            if( IS_CTCE_YAP( pDEVBLK->ctcexState ) )
            {

                // If available then we send the ATTN flag.
                if( IS_CTCE_YAV( pDEVBLK->ctcexState ) )
                   SET_CTCE_ATTN( pDEVBLK->ctcexState );

                // We always end up in the YWC state.
                SET_CTCE_YWC( pDEVBLK->ctcexState );

                // And then we always notify the other (y-)side.
                CTCT_Send( pDEVBLK, sCount, pIOBuf, pUnitStat, pResidual );

                // The above only proceeds when the matching SCB command
                // completes, unitistat and residual will be set,
                // and then we just return to available.
                SET_CTCE_YAV( pDEVBLK->ctcexState );
            }
            break;
        }
        // Non-Enhanced CTCT processing is retained.

        *pUnitStat = CSW_CE | CSW_DE;
        break;

    case 0x03:  // M0MMM011  NOP
        // -----------------------------------------------------------
        // CONTROL NO-OPERATON
        // -----------------------------------------------------------

        *pUnitStat = CSW_CE | CSW_DE;
        break;

    case 0x43:  // 00XXX011  SBM
        // -----------------------------------------------------------
        // SET BASIC MODE
        // -----------------------------------------------------------

        // Command reject if in basic mode
        if( pDEVBLK->ctcxmode == 0 )
        {
            pDEVBLK->sense[0] = SENSE_CR;
            *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;

            break;
        }

        // Reset extended mode and return normal status
        pDEVBLK->ctcxmode = 0;

        *pResidual = 0;
        *pUnitStat = CSW_CE | CSW_DE;

        break;

    case 0xC3:  // 11000011  SEM
        // -----------------------------------------------------------
        // SET EXTENDED MODE
        // -----------------------------------------------------------

        pDEVBLK->ctcxmode = 1;

        *pResidual = 0;
        *pUnitStat = CSW_CE | CSW_DE;

        break;

    case 0xE3:  // 11100011
        // -----------------------------------------------------------
        // PREPARE (PREP)
        // -----------------------------------------------------------

        *pUnitStat = CSW_CE | CSW_DE;

        break;

    case 0x14:  // XXX10100  SCB
        // -----------------------------------------------------------
        // SENSE COMMAND BYTE
        // -----------------------------------------------------------

        // Enhanced CTCT processing.
        if( pDEVBLK->ctcemode )
        {

            // If we were awaiting this SCB, i.e. if we are in the
            // YWC state, then we signal matching command received.
            if( IS_CTCE_YWC( pDEVBLK->ctcexState ) )
            {
                SET_CTCE_MATCH( pDEVBLK->ctcexState );

                // We complete the above signalling (if any) to the other (y-)side.
                CTCT_Send( pDEVBLK, sCount, NULL, pUnitStat, pResidual );

                // And only then to we return to state YAV.
                SET_CTCE_YAV( pDEVBLK->ctcexState );
            }

            *pIOBuf = pDEVBLK->ctceyCmdSCB;
            *pResidual = sCount - 1;
            *pUnitStat = CSW_CE | CSW_DE;

            if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
            {
                logmsg( _("HHCCT080I %4.4X: SCB executed: CB=%2.2X (x=%2.2X y=%2.2X)\n"),
                        pDEVBLK->devnum, *pIOBuf, pDEVBLK->ctcexState, pDEVBLK->ctceyState );
            }

            break;
        }
        // Non-Enhanced CTCT processing is retained.

        *pUnitStat = CSW_CE | CSW_DE;
        break;

    case 0x04:  // 00000100  SENSE
      // -----------------------------------------------------------
      // SENSE
      // -----------------------------------------------------------

        // Command reject if in basic mode
        if( pDEVBLK->ctcxmode == 0 )
        {
            pDEVBLK->sense[0] = SENSE_CR;
            *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
            break;
        }

        // Calculate residual byte count
        iNum = ( sCount < pDEVBLK->numsense ) ?
            sCount : pDEVBLK->numsense;

        *pResidual = sCount - iNum;

        if( sCount < pDEVBLK->numsense )
            *pMore = 1;

        // Copy device sense bytes to channel I/O buffer
        memcpy( pIOBuf, pDEVBLK->sense, iNum );

        // Clear the device sense bytes
        memset( pDEVBLK->sense, 0, sizeof( pDEVBLK->sense ) );

        // Return unit status
        *pUnitStat = CSW_CE | CSW_DE;

        break;

    case 0xE4:  //  11100100  SID
        // -----------------------------------------------------------
        // SENSE ID
        // -----------------------------------------------------------

        // Calculate residual byte count
        iNum = ( sCount < pDEVBLK->numdevid ) ?
            sCount : pDEVBLK->numdevid;

        *pResidual = sCount - iNum;

        if( sCount < pDEVBLK->numdevid )
            *pMore = 1;

        // Copy device identifier bytes to channel I/O buffer
        memcpy( pIOBuf, pDEVBLK->devid, iNum );

        // Return unit status
        *pUnitStat = CSW_CE | CSW_DE;

        break;

    default:
        // ------------------------------------------------------------
        // INVALID OPERATION
        // ------------------------------------------------------------

        // Set command reject sense byte, and unit check status
        pDEVBLK->sense[0] = SENSE_CR;
        *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
    }

    // Clear all flags when in state YAP
    if( IS_CTCE_YAP( pDEVBLK->ctcexState ) )
    {
        CLR_CTCE_ALLF( pDEVBLK->ctcexState );
    }

    release_lock( &pDEVBLK->lock );
}

// ====================================================================
// CTCT Support
// ====================================================================

//
// CTCT_Init
//

static int  CTCT_Init( DEVBLK *dev, int argc, char *argv[] )
{
    char           str[80];            // Thread name
    int            rc;                 // Return code
    int            mtu;                // MTU size (binary)
    int            lport;              // Listen port (binary)
    int            rport;              // Destination port (binary)
    char*          listenp;            // Listening port number
    char*          remotep;            // Destination port number
    char*          mtusize;            // MTU size (characters)
    char*          remaddr;            // Remote IP address
    struct in_addr ipaddr;             // Work area for IP address
    BYTE           c;                  // Character work area
    TID            tid;                // Thread ID for server
    TID            tid2;               // Thread ID for read thread
    u_int          ctceWrPort = 0;     // 0=read port, 1=write port
    int            ctceSmlBin;         // Small size (binary)
    char*          ctceSmlChr;         // Small size (characters)
    CTCG_PARMBLK   parm;               // Parameters for the server
    char           address[20]="";     // temp space for IP address

    dev->devtype = 0x3088;

    dev->ctctype = CTC_CTCT;

    SetSIDInfo( dev, 0x3088, 0x08, 0x3088, 0x01 );

    // Check for correct number of arguments
    if( (argc != 4) && (argc != 5) )
    {
        logmsg( _("HHCCT002E %4.4X: Incorrect number of parameters\n"),
                dev->devnum );
        return -1;
    }

    // The first argument is the listening port number
    listenp = *argv++;

    if( strlen( listenp ) > 5 ||
        sscanf( listenp, "%u%c", &lport, &c ) != 1 ||
        lport < 1024 || lport > 65534 )
    {
        logmsg( _("HHCCT003E %4.4X: Invalid port number: %s\n"),
                dev->devnum, listenp );
        return -1;
    }

    // The second argument is the IP address or hostname of the
    // remote side of the point-to-point link
    remaddr = *argv++;

    if( inet_aton( remaddr, &ipaddr ) == 0 )
    {
        struct hostent *hp;

        if( ( hp = gethostbyname( remaddr ) ) != NULL )
        {
            memcpy( &ipaddr, hp->h_addr, hp->h_length );
            strcpy( address, inet_ntoa( ipaddr ) );
            remaddr = address;
        }
        else
        {
            logmsg( _("HHCCT004E %4.4X: Invalid IP address %s\n"),
                    dev->devnum, remaddr );
            return -1;
        }
    }

    // The third argument is the destination port number
    remotep = *argv++;

    if( strlen( remotep ) > 5 ||
        sscanf( remotep, "%u%c", &rport, &c ) != 1 ||
        rport < 1024 || rport > 65534 )
    {
        logmsg( _("HHCCT005E %4.4X: Invalid port number: %s\n"),
                dev->devnum, remotep );
        return -1;
    }

    // The fourth argument is the maximum transmission unit (MTU) size
    mtusize = *argv;

    if( strlen( mtusize ) > 5 ||
        sscanf( mtusize, "%u%c", &mtu, &c ) != 1 ||
        mtu < 46 || mtu > 65536 )
    {
        logmsg( _("HHCCT006E %4.4X: Invalid MTU size %s\n"),
                dev->devnum, mtusize );
        return -1;
    }

    // Enhanced CTCT configuration requires additional checks.
    if( mtu >= CTCE_MTU_MIN )
    {
        // Enhanced CTCT mode is selected as bufsize >= CTCE_MTU_MIN (=32776).
        dev->ctcemode = 1;

        // Enhanced CTCT needs extended mode from the start.
        dev->ctcxmode = 1;

        // Enhanced CTCT needs even lport and rport numbers in the config.
        if( lport % 2 )
        {
            logmsg( _("HHCCT081E %4.4X: Enhanced CTCT (bufize >= 32776) lport number not even: %s\n"),
                    dev->devnum, listenp );
            return -1;
        }
        if( rport % 2 )
        {
            logmsg( _("HHCCT082E %4.4X: Enhanced CTCT (bufize >= 32776) rport number not even: %s\n"),
                    dev->devnum, remotep );
            return -1;
        }

        // Mark both socket file descriptors as not yet connected.
        dev->fd = -1;
        dev->ctcefd = -1;

        // Enhanced CTCT only supports an optional 5th parameter,
        // the Small MTU size, which defaults to the MTU size.
        ctceSmlBin = mtu;
        if( argc == 5 )
        {
            ctceSmlChr = *(++argv);

            if( strlen( ctceSmlChr ) > 5 ||
                sscanf( ctceSmlChr, "%u%c", &ctceSmlBin, &c ) != 1 ||
                ctceSmlBin < 46 || ctceSmlBin > mtu )
            {
                ctceSmlBin = mtu;
                logmsg( _("HHCCT083W %4.4X: Invalid Small MTU size %s ignored\n"),
                        dev->devnum, ctceSmlChr );
            }
        }
        dev->ctceSndSml = ctceSmlBin;

    }

    // Set the device buffer size equal to the MTU size
    dev->bufsize = mtu;

    // Initialize the file descriptor for the socket connection

    // It's a little confusing, but we're using a couple of the
    // members of the server paramter structure to initiate the
    // outgoing connection.  Saves a couple of variable declarations,
    // though.  If we feel strongly about it, we can declare separate
    // variables...

    // Enhanced CTCT will require a second pass for the odd port number.
    for(ctceWrPort = 0; ctceWrPort <= dev->ctcemode; ctceWrPort++ )
    {
        // make a TCP socket
        parm.listenfd[ctceWrPort] = socket( AF_INET, SOCK_STREAM, 0 );

        if( parm.listenfd[ctceWrPort] < 0 )
        {
            logmsg( _("HHCCT007E %4.4X: Error creating socket: %s\n"),
                    dev->devnum, strerror( HSO_errno ) );
            CTCX_Close( dev );
            return -1;
        }

        // bind socket to our local port
        // (might seem like overkill, and usually isn't done, but doing this
        // bind() to the local port we configure gives the other end a chance
        // at validating the connection request)
        memset( &(parm.addr), 0, sizeof( parm.addr ) );
        parm.addr.sin_family      = AF_INET;
        parm.addr.sin_port        = htons(lport + ctceWrPort);
        parm.addr.sin_addr.s_addr = htonl(INADDR_ANY);

        rc = bind( parm.listenfd[ctceWrPort],
                   (struct sockaddr *)&parm.addr,
                   sizeof( parm.addr ) );
        if( rc < 0 )
        {
            logmsg( _("HHCCT008E %4.4X: Error binding to socket (port %d): %s\n"),
                    dev->devnum, lport + ctceWrPort, strerror( HSO_errno ) );
            CTCX_Close( dev );
            return -1;
        }

        // initiate a connection to the other end
        memset( &(parm.addr), 0, sizeof( parm.addr ) );
        parm.addr.sin_family = AF_INET;
        // the even (=read) port must connect to the odd (=write) port
        // at the other side and vice-versa
        parm.addr.sin_port   = htons(rport + (ctceWrPort + 1)%2);
        parm.addr.sin_addr   = ipaddr;
        rc = connect( parm.listenfd[ctceWrPort],
                      (struct sockaddr *)&parm.addr,
                      sizeof( parm.addr ) );

        // if connection was not successful, start a server
        if( rc < 0 )
        {
            // used to pass parameters to the server thread
            CTCG_PARMBLK* arg;

            logmsg( _("HHCCT009I %4.4X: Connect from %d to %s:%d failed, starting server\n"),
                    dev->devnum, lport + ctceWrPort, remaddr, rport + (ctceWrPort + 1)%2 );

            // probably don't need to do this, not sure...
            close_socket( parm.listenfd[ctceWrPort] );

            parm.listenfd[ctceWrPort] = socket( AF_INET, SOCK_STREAM, 0 );

            if( parm.listenfd[ctceWrPort] < 0 )
            {
                logmsg( _("HHCCT010E %4.4X: Error creating socket: %s\n"),
                        dev->devnum, strerror( HSO_errno ) );
                CTCX_Close( dev );
                return -1;
            }

            // set up the listening port
            memset( &(parm.addr), 0, sizeof( parm.addr ) );

            parm.addr.sin_family      = AF_INET;
            parm.addr.sin_port        = htons(lport + ctceWrPort) ;
            parm.addr.sin_addr.s_addr = htonl(INADDR_ANY);

            if( bind( parm.listenfd[ctceWrPort],
                      (struct sockaddr *)&parm.addr,
                      sizeof( parm.addr ) ) < 0 )
            {
                logmsg( _("HHCCT011E %4.4X: Error binding to socket (port=%d): %s\n"),
                        dev->devnum, lport + ctceWrPort, strerror( HSO_errno ) );
                CTCX_Close( dev );
                return -1;
            }

            if( listen( parm.listenfd[ctceWrPort], 1 ) < 0 )
            {
                logmsg( _("HHCCT012E %4.4X: Error on call to listen (port=%d): %s\n"),
                        dev->devnum, lport + ctceWrPort, strerror( HSO_errno ) );
                CTCX_Close( dev );
                return -1;
            }

            // we are listening, so create a thread to accept connection
            arg = malloc( sizeof( CTCG_PARMBLK ) );
            memcpy( arg, &parm, sizeof( parm ) );
            arg->dev = dev;
            arg->ctceWrPort = ctceWrPort;
            snprintf(str,sizeof(str),"CTCT %4.4X ListenThread %d",dev->devnum, ctceWrPort);
            str[sizeof(str)-1]=0;
            create_thread( &tid, JOINABLE, CTCT_ListenThread, arg, str );

        }
        else  // successfully connected (outbound) to the other end
        {
            logmsg( _("HHCCT013I %4.4X: Connect port %d to %s:%d succes\n"),
                    dev->devnum, lport + ctceWrPort, remaddr, rport + (ctceWrPort + 1)%2 );

            // The even local port (form the config) is for writing
            if( ctceWrPort == 0 )
            {
                dev->fd = parm.listenfd[ctceWrPort];
            }
            else
            {
                // The next odd local port (form the config) is for reading
                dev->ctcefd = parm.listenfd[ctceWrPort];

                // This side is ready to start receiving and sending so we
                // start a read thread to do the receiving part; identical
                // code will be found in the CTCT_ListenThread after a
                // successful connect was accepted there.
                snprintf(str,sizeof(str),"CTCT %4.4X RecvThread %d",dev->devnum, ctceWrPort);
                str[sizeof(str)-1]=0;
                create_thread( &tid2, JOINABLE, CTCT_RecvThread, dev, str );
            }
        }

        // for cosmetics, since we are successfully connected or serving,
        // fill in some details for the panel.
        // Also used for connection verification in CTCT_ListenThread
        sprintf( dev->filename, "%s:%d", remaddr, rport + (ctceWrPort + 1)%2 );
    }

    // Enhanced CTCT adapter intiialization for command register and CB.
    dev->ctcexCmd = 0x00;
    dev->ctceyCmd = 0x00;
    dev->ctceyCmdSCB = 0x00;

    // Enhanced CTCT adapter sides are state-aware, with initial
    // state "Available" = YAV which corresponds to column 5 in
    // the table 2.13 in SA22-7203-00, i.e. we consider both
    // x- and y-side READY from the start. ALL Flags are cleared.
    CLR_CTCE_ALLF( dev->ctcexState );
//  SET_CTCE_NRDY( dev->ctcexState ); // Please see also earlier notes.
    SET_CTCE_YAV ( dev->ctcexState );
    CLR_CTCE_ALLF( dev->ctceyState );
//  SET_CTCE_NRDY( dev->ctceyState ); // Please see also earlier notes.
    SET_CTCE_YAV ( dev->ctceyState );

    // Initialize the 12 bits Send->Recv packet counter with
    // the leftmost 4 bits C of the CCUU devnum at this side,
    // which helps distinguishing same-host traffic if the
    // Send-Recv side CCUU's have a diffent leftmost C.
    dev->ctcePktCnt = dev->devnum & 0xF000;

    // Initialize the CTC lock and condition used to signal
    // reception of a command matching the dependent one.
    initialize_lock( &dev->ctceEventLock );
    initialize_condition( &dev->ctceEvent );

    return 0;
}

//
// CTCT_Send
//

static void   CTCT_Send( DEVBLK* pDEVBLK,   U16   sCount,
                         BYTE*   pIOBuf,    BYTE* pUnitStat,
                         U16*    pResidual )
{
    CTCE_SOKPFX   *pSokBuf;                 // overlay for buf inside DEVBLK
    int            rc;                      // Return code

    int            i;                       // temp counter
    U32            XORChk = 0;              // XOR of sent buffer for checking
    BYTE          *pXOR = (BYTE*)&XORChk;   // -> XORChk
    BYTE          *pBuf = pDEVBLK->buf;     // temp pointer inside buf

    pSokBuf = (CTCE_SOKPFX*) pDEVBLK->buf;
    pSokBuf->CmdReg = pDEVBLK->ctcexCmd;
    pSokBuf->FsmSta = pDEVBLK->ctcexState;
    pSokBuf->sCount = sCount;
    pSokBuf->PktCnt = ++pDEVBLK->ctcePktCnt;
    pSokBuf->SndLen = pDEVBLK->ctceSndSml;

    // We only ever Send if the sockets are connected.
    if( ( pDEVBLK->fd < 0) || ( pDEVBLK->ctcefd < 0) )
        return ;

    // Only a (non-WEOF) write command data includes sending the IOBuf.
//  if( IS_CCW_WRITE_NOEOF( pDEVBLK->ctcexCmd ) )
    if( IS_CTCE_CCW_WRT( pDEVBLK->ctcexCmd ) )
    {
        memcpy( pDEVBLK->buf + sizeof(CTCE_SOKPFX), pIOBuf, sCount );

        // Increase the SndLen if the sCount is too large.
        if( pSokBuf->SndLen < ( sCount + sizeof(CTCE_SOKPFX) ) )
            pSokBuf->SndLen = ( sCount + sizeof(CTCE_SOKPFX) );


//1     if( pDEVBLK->ccwstep )
//1         packet_trace( pIOBuf, sCount );

        // If bufsize (init from the MTU parameter) is not large enough
        // then we will have a severe error as the CTC will not connect.
        if( pDEVBLK->bufsize < pSokBuf->SndLen )
        {
            logmsg( _("HHCCT084S %4.4X: bufsize parameter %d is too small; increase at least to %d\n"),
                    pDEVBLK->devnum, pDEVBLK->bufsize, pSokBuf->SndLen );
        }
    }

    // Write the all of this to the other (y-)side.
    rc = write_socket( pDEVBLK->fd, pDEVBLK->buf, pSokBuf->SndLen );

    // Trace the IP packet just sent if needed.
    if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
    {
        XORChk = 0;
        pBuf = pDEVBLK->buf;
        for(i = 0; i < pSokBuf->SndLen; i++)
        {
           if( (i % 4) == 0 )
               pXOR = (BYTE*)&XORChk;
           *pXOR++ ^= *pBuf++;
        }
        logmsg( _("HHCCT085I %4.4X: Send %4.4X->%s %s=%2.2X x=%2.2X y=%2.2X l=%4.4X k=%8.8X\n"),
                pDEVBLK->devnum, pSokBuf->PktCnt, pDEVBLK->filename,
                CTCE_CmdStr[CTCE_Cmd[pDEVBLK->ctcexCmd]], pDEVBLK->ctcexCmd,
                pDEVBLK->ctcexState, pDEVBLK->ctceyState, sCount, XORChk ) ;
        if( pDEVBLK->ccwtrace )
            packet_trace( pDEVBLK->buf, pSokBuf->SndLen );
    }

    if( rc < 0 )
    {
        logmsg( _("HHCCT086E %4.4X: Error writing to %s: %s\n"),
                pDEVBLK->devnum, pDEVBLK->filename,
                strerror( HSO_errno ) );

        pDEVBLK->sense[0] = SENSE_EC;
        *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
        return;
    }

    // The ATTN flag signals entering a Working(D) state.  If this was
    // because of a READ command then we need to wait for the matching
    // WRITE command from the other (y-)side to arrive.

    // Please note that initially we believed such processing to be
    // needed for READ and CONTROL commands as well, but this turned out
    // not to be the case.  Even stronger, we have not yet seen a READ
    // command BEFORE the matching WRITE command, thus this processing
    // is perhaps never needed at all ...
    //                                                 (PJJ, April 2014)
    else if( IS_CTCE_ATTN( pDEVBLK->ctcexState ) &&
             IS_CCW_READ( pDEVBLK->ctcexCmd ) )
    {
        struct timespec waittime;
        struct timeval  now;

        // Any request to signal attention has now been sent.
        CLR_CTCE_ATTN( pDEVBLK->ctcexState );

        gettimeofday( &now, NULL );

        waittime.tv_sec  = now.tv_sec  + CTC_READ_TIMEOUT_SECS;
        waittime.tv_nsec = now.tv_usec * 1000;

        obtain_lock( &pDEVBLK->ctceEventLock );

        rc = timed_wait_condition( &pDEVBLK->ctceEvent,
                                   &pDEVBLK->ctceEventLock,
                                   &waittime );
        release_lock( &pDEVBLK->ctceEventLock );

        // Trace the reception of the matching command, or the wait timeout (RC=138).
        if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
            logmsg( _("HHCCT087I %4.4X: Send %4.4X to %s %s=%2.2X x=%2.2X y=%2.2X: wait RC=%d\n"),
                    pDEVBLK->devnum, pSokBuf->PktCnt, pDEVBLK->filename,
                    CTCE_CmdStr[CTCE_Cmd[pDEVBLK->ctcexCmd]], pDEVBLK->ctcexCmd,
                    pDEVBLK->ctcexState, pDEVBLK->ctceyState, rc ) ;

        // First we check for Halt or Clear Subchannel
        if( rc == ETIMEDOUT || rc == EINTR )
        {
            // check for halt condition
            if( pDEVBLK->scsw.flag2 & SCSW2_FC_HALT ||
                pDEVBLK->scsw.flag2 & SCSW2_FC_CLEAR )
            {
                if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
                    logmsg( _("HHCCT088I %4.4X: Halt or Clear Recognized\n"),
                            pDEVBLK->devnum );

                *pUnitStat = CSW_CE | CSW_DE;
                *pResidual = sCount;
            }

            // Other timeouts or errors should not occur.
            else
            {
                *pUnitStat = CSW_CE | CSW_DE | CSW_UC | CSW_SM;
                pDEVBLK->sense[0] = 0;
            }
            return;
        }

        // A WRITE EOF command from the other side will have resulted
        // in the WEOF flag being set.  If this was a matching command
        // for a READ then unit exception needs to be included.
        else if( IS_CTCE_WEOF( pDEVBLK->ctcexState ) )
        {
            *pResidual = 0;
            *pUnitStat  = CSW_CE | CSW_DE | CSW_UX;

            if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
                logmsg( _("HHCCT089I %4.4X: Recv %4.4X to %s %s=%2.2X x=%2.2X y=%2.2X: WEOF -> UX\n"),
                        pDEVBLK->devnum, pSokBuf->PktCnt, pDEVBLK->filename,
                        CTCE_CmdStr[CTCE_Cmd[pDEVBLK->ctcexCmd]], pDEVBLK->ctcexCmd,
                        pDEVBLK->ctcexState, pDEVBLK->ctceyState ) ;
            return;
        }
    }

    // Reset the attention signal flag.
    CLR_CTCE_ATTN( pDEVBLK->ctcexState );

    // If the command (by now matched) was a READ command, then the
    // other (y-)side data is available in the DEVBLK buf, so we
    // can copy it into the IO channel buffer and compute residual.
    if( IS_CCW_READ( pDEVBLK->ctcexCmd ) )
    {
        memcpy( pIOBuf, pDEVBLK->buf + sizeof(CTCE_SOKPFX),
            ( pSokBuf->sCount <= pDEVBLK->bufsize) ? pSokBuf->sCount : pDEVBLK->bufsize );
        *pResidual = sCount - pSokBuf->sCount;

//1     if( pDEVBLK->ccwstep )
//1         packet_trace( pIOBuf, sCount );
    }
    else
        *pResidual = 0;

    if( IS_CCW_WRITE( pDEVBLK->ctcexCmd ) )
        *pUnitStat = 0;
    else
        *pUnitStat = CSW_CE | CSW_DE;
    return;
}

//
// CTCT_Write
//

static void  CTCT_Write( DEVBLK* pDEVBLK,   U16   sCount,
                         BYTE*   pIOBuf,    BYTE* pUnitStat,
                         U16*    pResidual )
{
    PCTCIHDR   pFrame;                  // -> Frame header
    PCTCISEG   pSegment;                // -> Segment in buffer
    U16        sOffset;                 // Offset of next frame
    U16        sSegLen;                 // Current segment length
    U16        sDataLen;                // Length of IP Frame data
    int        iPos;                    // Offset into buffer
    U16        i;                       // Array subscript
    int        rc;                      // Return code
    BYTE       szStackID[33];           // VSE IP stack identity
    U32        iStackCmd;               // VSE IP stack command

    // Check that CCW count is sufficient to contain block header
    if( sCount < sizeof( CTCIHDR ) )
    {
        logmsg( _("HHCCT014E %4.4X: Write CCW count %u is invalid\n"),
                pDEVBLK->devnum, sCount );

        pDEVBLK->sense[0] = SENSE_DC;
        *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
        return;
    }

    // Fix-up frame pointer
    pFrame = (PCTCIHDR)pIOBuf;

    // Extract the frame length from the header
    FETCH_HW( sOffset, pFrame->hwOffset );


    // Check for special VSE TCP/IP stack command packet
    if( sOffset == 0 && sCount == 40 )
    {
        // Extract the 32-byte stack identity string
        for( i = 0;
             i < sizeof( szStackID ) - 1 && i < sCount - 4;
             i++)
            szStackID[i] = guest_to_host( pIOBuf[i+4] );

        szStackID[i] = '\0';

        // Extract the stack command word
        FETCH_FW( iStackCmd, *((FWORD*)&pIOBuf[36]) );

        // Display stack command and discard the packet
        logmsg( _("HHCCT015I %4.4X: Interface command: %s %8.8X\n"),
                pDEVBLK->devnum, szStackID, iStackCmd );

        *pUnitStat = CSW_CE | CSW_DE;
        *pResidual = 0;
        return;
    }

    // Check for special L/390 initialization packet
    if( sOffset == 0 )
    {
        // Return normal status and discard the packet
        *pUnitStat = CSW_CE | CSW_DE;
        *pResidual = 0;
        return;
    }

#if 0
    // Notes: It appears that TurboLinux has gotten sloppy in their
    //        ways. They are now giving us buffer sizes that are
    //        greater than the CCW count, but the segment size
    //        is within the count.
    // Check that the frame offset is valid
    if( sOffset < sizeof( CTCIHDR ) || sOffset > sCount )
    {
        logmsg( _("CTC101W %4.4X: Write buffer contains invalid "
                  "frame offset %u\n"),
                pDEVBLK->devnum, sOffset );

        pDEVBLK->sense[0] = SENSE_CR;
        *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
        return;
    }
#endif

    // Adjust the residual byte count
    *pResidual -= sizeof( CTCIHDR );

    // Process each segment in the buffer
    for( iPos  = sizeof( CTCIHDR );
         iPos  < sOffset;
         iPos += sSegLen )
    {
        // Check that the segment is fully contained within the block
        if( iPos + sizeof( CTCISEG ) > sOffset )
        {
            logmsg( _("HHCCT016E %4.4X: Write buffer contains incomplete "
                      "segment header at offset %4.4X\n"),
                    pDEVBLK->devnum, iPos );

            pDEVBLK->sense[0] = SENSE_DC;
            *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
            return;
        }

        // Fix-up segment header in the I/O buffer
        pSegment = (PCTCISEG)(pIOBuf + iPos);

        // Extract the segment length from the segment header
        FETCH_HW( sSegLen, pSegment->hwLength );

        // Check that the segment length is valid
        if( ( sSegLen        < sizeof( CTCISEG ) ) ||
            ( iPos + sSegLen > sOffset           ) ||
            ( iPos + sSegLen > sCount            ) )
        {
            logmsg( _("HHCCT017E %4.4X: Write buffer contains invalid "
                    "segment length %u at offset %4.4X\n"),
                    pDEVBLK->devnum, sSegLen, iPos );

            pDEVBLK->sense[0] = SENSE_DC;
            *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
            return;
        }

        // Calculate length of IP frame data
        sDataLen = sSegLen - sizeof( CTCISEG );

        // Trace the IP packet before sending
        if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
        {
            logmsg( _("HHCCT018I %4.4X: Sending packet to %s:\n"),
                    pDEVBLK->devnum, pDEVBLK->filename );
            if( pDEVBLK->ccwtrace )
                packet_trace( pSegment->bData, sDataLen );
        }

        // Write the IP packet
        rc = write_socket( pDEVBLK->fd, pSegment->bData, sDataLen );

        if( rc < 0 )
        {
            logmsg( _("HHCCT019E %4.4X: Error writing to %s: %s\n"),
                    pDEVBLK->devnum, pDEVBLK->filename,
                    strerror( HSO_errno ) );

            pDEVBLK->sense[0] = SENSE_EC;
            *pUnitStat        = CSW_CE | CSW_DE | CSW_UC;
            return;
        }

        // Adjust the residual byte count
        *pResidual -= sSegLen;

        // We are done if current segment satisfies CCW count
        if( iPos + sSegLen == sCount )
        {
            *pResidual -= sSegLen;
            *pUnitStat = CSW_CE | CSW_DE;
            return;
        }
    }

    // Set unit status and residual byte count
    *pUnitStat = CSW_CE | CSW_DE;
    *pResidual = 0;
}

//
// CTCT_Read
//

static void  CTCT_Read( DEVBLK* pDEVBLK,   U16   sCount,
                        BYTE*   pIOBuf,    BYTE* pUnitStat,
                        U16*    pResidual, BYTE* pMore )
{
    PCTCIHDR    pFrame   = NULL;       // -> Frame header
    PCTCISEG    pSegment = NULL;       // -> Segment in buffer
    fd_set      rfds;                  // Read FD_SET
    int         iRetVal;               // Return code from 'select'
    ssize_t     iLength  = 0;

    static struct timeval tv;          // Timeout time for 'select'


    // Limit how long we should wait for data to come in
    FD_ZERO( &rfds );
    FD_SET( pDEVBLK->fd, &rfds );

    tv.tv_sec  = CTC_READ_TIMEOUT_SECS;
    tv.tv_usec = 0;

    iRetVal = select( pDEVBLK->fd + 1, &rfds, NULL, NULL, &tv );

    switch( iRetVal )
    {
    case 0:
        *pUnitStat = CSW_CE | CSW_DE | CSW_UC | CSW_SM;
        pDEVBLK->sense[0] = 0;
        return;

    case -1:
        if( HSO_errno == HSO_EINTR )
            return;

        logmsg( _("HHCCT020E %4.4X: Error reading from %s: %s\n"),
                pDEVBLK->devnum, pDEVBLK->filename, strerror( HSO_errno ) );

        pDEVBLK->sense[0] = SENSE_EC;
        *pUnitStat = CSW_CE | CSW_DE | CSW_UC;
        return;

    default:
        break;
    }

    // Read an IP packet from the TUN device
    iLength = read_socket( pDEVBLK->fd, pDEVBLK->buf, pDEVBLK->bufsize );

    // Check for other error condition
    if( iLength < 0 )
    {
        logmsg( _("HHCCT021E %4.4X: Error reading from %s: %s\n"),
                pDEVBLK->devnum, pDEVBLK->filename, strerror( HSO_errno ) );
        pDEVBLK->sense[0] = SENSE_EC;
        *pUnitStat = CSW_CE | CSW_DE | CSW_UC;
        return;
    }

    // Trace the packet received from the TUN device
    if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
    {
        logmsg( _("HHCCT022I %4.4X: Received packet from %s (%d bytes):\n"),
                pDEVBLK->devnum, pDEVBLK->filename, iLength );
        packet_trace( pDEVBLK->buf, iLength );
    }

    // Fix-up Frame pointer
    pFrame = (PCTCIHDR)pIOBuf;

    // Fix-up Segment pointer
    pSegment = (PCTCISEG)( pIOBuf + sizeof( CTCIHDR ) );

    // Initialize segment
    memset( pSegment, 0, iLength + sizeof( CTCISEG ) );

    // Update next frame offset
    STORE_HW( pFrame->hwOffset,
              iLength + sizeof( CTCIHDR ) + sizeof( CTCISEG ) );

    // Store segment length
    STORE_HW( pSegment->hwLength, iLength + sizeof( CTCISEG ) );

    // Store Frame type
    STORE_HW( pSegment->hwType, ETH_TYPE_IP );

    // Copy data
    memcpy( pSegment->bData, pDEVBLK->buf, iLength );

    // Fix-up frame pointer and terminate block
    pFrame = (PCTCIHDR)( pIOBuf + sizeof( CTCIHDR ) +
                         sizeof( CTCISEG ) + iLength );
    STORE_HW( pFrame->hwOffset, 0x0000 );

    // Calculate #of bytes returned including two slack bytes
    iLength += sizeof( CTCIHDR ) + sizeof( CTCISEG ) + 2;

    if( sCount < iLength )
    {
        *pMore     = 1;
        *pResidual = 0;

        iLength    = sCount;
    }
    else
    {
        *pMore      = 0;
        *pResidual -= iLength;
    }

    // Set unit status
    *pUnitStat = CSW_CE | CSW_DE;
}

//
// CTCT_ListenThread
//

static void*  CTCT_ListenThread( void* argp )
{
    int          connfd;
    socklen_t    servlen;
    char         str[80];
    CTCG_PARMBLK parm;
    TID          tid2;               // Thread ID for read thread

    // set up the parameters passed via create_thread
    parm = *((CTCG_PARMBLK*) argp);
    free( argp );

    for( ; ; )
    {
        servlen = sizeof(parm.addr);

        // await a connection
        connfd = accept( parm.listenfd[parm.ctceWrPort],
                         (struct sockaddr *)&parm.addr,
                         &servlen );

        sprintf( str, "%s:%d",
                 inet_ntoa( parm.addr.sin_addr ),
                 ntohs( parm.addr.sin_port ) - (parm.ctceWrPort + 1)%2 );

        if( strcmp( str, parm.dev->filename ) != 0 )
        {
            logmsg( _("HHCCT023E %4.4X: Incorrect client or config error\n"
                      "                Config=%s+%d, connecting client=%s\n"),
                    parm.dev->devnum,
                    parm.dev->filename, parm.ctceWrPort, str);
            close_socket( connfd );
        }
        else
        {
            // The even local port (as in the config) is for writing
            if( parm.ctceWrPort == 0 )
            {
                parm.dev->fd = connfd;
            }
            else
            {
                // The next odd local port is for reading
                parm.dev->ctcefd = connfd;

                // This side is ready to start receiving and sending so we
                // start a read thread to do the receiving part; identical
                // code will be found in the CTCT_ListenThread after a
                // successful connect was accepted there
                snprintf(str,sizeof(str),"CTCT %4.4X RecvThread %d",parm.dev->devnum, parm.ctceWrPort);
                str[sizeof(str)-1]=0;
                create_thread( &tid2, JOINABLE, CTCT_RecvThread, parm.dev, str );
            }
        }

        // Ok, so having done that we're going to loop back to the
        // accept().  This was meant to handle the connection failing
        // at the other end; this end will be ready to accept another
        // connection.  Although this will happen, I'm sure you can
        // see the possibility for bad things to occur (eg if another
        // Hercules tries to connect).  This will also be fixed RSN.
    }

    return NULL;    // make compiler happy
}

//
// CTCT_RecvThread
//

static void*  CTCT_RecvThread( DEVBLK* pDEVBLK )
{
    CTCE_SOKPFX   *pSokBuf;                      // overlay for buf inside DEVBLK
    ssize_t        iLength  = 0;
    BYTE          *buf;                          //-> Device recv data buffer

    U64            ctcePktCnt = 0;               // Recvd Packet Count
    U64            ctceBytCnt = 0;               // Recvd Byte Count

    int            i;                            // temp counter
    U32            XORChk = 0;                   // XOR of sent buffer for checking
    BYTE          *pXOR = (BYTE*)&XORChk;        // -> XORChk
    BYTE          *pBuf = pDEVBLK->buf;          // temp pointer inside buf

    logmsg( _("HHCCT090I %4.4X: Started %s reader thread (bufsize=%d,%d)\n"),
            pDEVBLK->devnum, pDEVBLK->filename, pDEVBLK->bufsize, pDEVBLK->ctceSndSml );

    // avoid having to lock the DEVBLK whilst awaiting data to arrive via read_socket
    buf = malloc( pDEVBLK->bufsize );
    pSokBuf = (CTCE_SOKPFX*)buf;

    // This thread will loop until we receive a zero-length packet
    for( ; ; )
    {
        // We read whatever the other (y-)side of the CTC has sent us,
        // which by now won't block until the complete bufsize is received.
//      iLength = read_socket( pDEVBLK->ctcefd, buf, CTCE_SND_SML );
        iLength = read_socket( pDEVBLK->ctcefd, buf, pDEVBLK->ctceSndSml );

        // Followed by the receiving the rest if the default SndLen was too small.
        if( ( pDEVBLK->ctceSndSml < pSokBuf->SndLen ) && ( iLength != 0 ) )
            iLength += read_socket( pDEVBLK->ctcefd, buf + pDEVBLK->ctceSndSml,
                pSokBuf->SndLen - pDEVBLK->ctceSndSml );

//      if( CTCE_SND_SML < ((CTCE_SOKPFX*)buf)->SndLen && iLength != 0 )
//          iLength += read_socket( pDEVBLK->ctcefd, buf + CTCE_SND_SML,
//              ((CTCE_SOKPFX*)buf)->SndLen - CTCE_SND_SML );

        // In case we are closing down this thread can end.
        if( iLength == 0 )
        {
            logmsg( _("HHCCT091E %4.4X: Zero length read from %s\n"),
                    pDEVBLK->devnum, pDEVBLK->filename );
            logmsg( _("HHCCT092E %4.4X: %d MB received in %d packets\n"),
                    pDEVBLK->devnum, ctceBytCnt / 1048576 , ctcePktCnt );
            free( buf );
            return NULL;    // make compiler happy
        }

        // Changes to DEVBLK must be lock protected as other threads might update as well.
        obtain_lock( &pDEVBLK->lock );

        // Check for other error condition
        if( iLength < 0 )
        {
            logmsg( _("HHCCT093E %4.4X: Error reading from %s: %s\n"),
                    pDEVBLK->devnum, pDEVBLK->filename, strerror( HSO_errno ) );
            pDEVBLK->sense[0] = SENSE_EC;
            pDEVBLK->scsw.unitstat = CSW_CE | CSW_DE | CSW_UC;
        }
        else
        {

            ctcePktCnt += 1 ;
            ctceBytCnt += iLength ;
//          pSokBuf = (CTCE_SOKPFX*)buf;
            pDEVBLK->ctceyCmd =  pSokBuf->CmdReg;
            pDEVBLK->ctceyState = pSokBuf->FsmSta;

            // Trace the packet received from the other side of the CTC
            if( pDEVBLK->ccwtrace || pDEVBLK->ccwstep )
            {
                XORChk = 0;
                pBuf = buf;
                for(i = 0; i < iLength; i++)
                {
                   if( (i % 4) == 0 )
                       pXOR = (BYTE*)&XORChk;
                   *pXOR++ ^= *pBuf++;
                }
                logmsg( _("HHCCT094I %4.4X: Recv %4.4X<-%s %s=%2.2X x=%2.2X y=%2.2X l=%4.4X k=%8.8X\n"),
                        pDEVBLK->devnum, pSokBuf->PktCnt, pDEVBLK->filename,
                        CTCE_CmdStr[CTCE_Cmd[pDEVBLK->ctceyCmd]], pDEVBLK->ctceyCmd,
                        pDEVBLK->ctcexState, pDEVBLK->ctceyState, pSokBuf->sCount, XORChk ) ;
                if( pDEVBLK->ccwtrace )
                    packet_trace( buf, iLength );
            }

            // Only if the other (y-)side sent us a write command will
            // we copy the socket buffer into the device buffer.
            if( IS_CTCE_CCW_WRT( pDEVBLK->ctceyCmd ) )
                memcpy( pDEVBLK->buf + sizeof(CTCE_SOKPFX), buf + sizeof(CTCE_SOKPFX), pSokBuf->sCount );

            // If the other side sent us a WRITE EOF command
            // then we just set the WEOF flag on our side.
            if( IS_CTCE_CCW_WEF(  pDEVBLK->ctceyCmd ) )
            {
                SET_CTCE_WEOF( pDEVBLK->ctcexState );

                // But only if this does NOT match a Working(D) Read command will it remain set.
                if( ( pDEVBLK->ccwtrace || pDEVBLK->ccwstep ) && (!IS_CTCE_MATCH( pDEVBLK->ctceyState ) ) )
                    logmsg( _("HHCCT095I %4.4X: Recv %4.4X<-%s %s=%2.2X x=%2.2X y=%2.2X: WEOF ->set\n"),
                            pDEVBLK->devnum, pSokBuf->PktCnt, pDEVBLK->filename,
                            CTCE_CmdStr[CTCE_Cmd[pDEVBLK->ctceyCmd]], pDEVBLK->ctceyCmd,
                            pDEVBLK->ctcexState, pDEVBLK->ctceyState ) ;
            }

            // If we were sent this command because the other side was
            // not ready yet, then we need to signal a device end and
            // reset SENSE[0] bits 1 (SENSE_IR) and 7 (SENSE_OC).
            // However, as stated earlier, this NOT READY handling did
            // not work and was effectively disabled by initializing
            // both adapter sides as being READY from the start.
            if( IS_CTCE_NRDY( pDEVBLK->ctceyState ) )
            {
                pDEVBLK->sense[0] &= ~(SENSE_IR | SENSE_OC);
                device_attention( pDEVBLK, CSW_DE );
            }

            // If the other (y-)side sent us a command that put them
            // in the Working(D) state asking us to signal attention
            // then we copy their state into ours.  We also copy the
            // y command register for the SCB command to read from.
            else if( IS_CTCE_ATTN( pDEVBLK->ctceyState ) )
            {
                CLR_CTCE_ATTN( pDEVBLK->ctceyState );
                pDEVBLK->ctcexState = pDEVBLK->ctceyState;
                pDEVBLK->ctceyCmdSCB = pDEVBLK->ctceyCmd;

                // and then we signal that attention.
                device_attention( pDEVBLK, CSW_ATTN );
            }

            // If the other (y-)side sent us a matching command then
            // the Working(D) state will cease into YAV state but
            // CTCT_Send is awaiting us to signal that condition.
            // We clear the SCB command buffer input.
            else if( IS_CTCE_MATCH( pDEVBLK->ctceyState ) )
            {
                CLR_CTCE_MATCH( pDEVBLK->ctceyState );
                pDEVBLK->ctceyCmdSCB = 0;
                obtain_lock( &pDEVBLK->ctceEventLock );
                signal_condition( &pDEVBLK->ctceEvent );
                release_lock( &pDEVBLK->ctceEventLock );
            }

        }

        release_lock( &pDEVBLK->lock );
    }

    free( buf );

    return NULL;    // make compiler happy
}

// ====================================================================
// VMNET Support -- written by Willem Konynenberg
// ====================================================================

/*-------------------------------------------------------------------*/
/* Definitions for SLIP encapsulation                                */
/*-------------------------------------------------------------------*/
#define SLIP_END        0300
#define SLIP_ESC        0333
#define SLIP_ESC_END    0334
#define SLIP_ESC_ESC    0335

/*-------------------------------------------------------------------*/
/* Functions to support vmnet written by Willem Konynenberg          */
/*-------------------------------------------------------------------*/
static int start_vmnet(DEVBLK *dev, DEVBLK *xdev, int argc, char *argv[])
{
int sockfd[2];
int r, i;
char *ipaddress;

    if (argc < 2) {
        logmsg (_("HHCCT024E %4.4X: Not enough arguments to start vmnet\n"),
                        dev->devnum);
        return -1;
    }

    ipaddress = argv[0];
    argc--;
    argv++;

    if (socketpair (AF_UNIX, SOCK_STREAM, 0, sockfd) < 0) {
        logmsg (_("HHCCT025E %4.4X: Failed: socketpair: %s\n"),
                        dev->devnum, strerror(errno));
        return -1;
    }

    r = fork ();

    if (r < 0) {
        logmsg (_("HHCCT026E %4.4X: Failed: fork: %s\n"),
                        dev->devnum, strerror(errno));
        return -1;
    } else if (r == 0) {
        /* child */
        close (0);
        close (1);
        dup (sockfd[1]);
        dup (sockfd[1]);
        r = (sockfd[0] > sockfd[1]) ? sockfd[0] : sockfd[1];
        for (i = 3; i <= r; i++) {
            close (i);
        }

        /* the ugly cast is to silence a compiler warning due to const */
        execv (argv[0], (EXECV_ARG2_ARGV_T)argv);

        exit (1);
    }

    close (sockfd[1]);
    dev->fd = sockfd[0];
    xdev->fd = sockfd[0];

    /* We just blindly copy these out in the hope vmnet will pick them
     * up correctly.  I don't feel like implementing a complete login
     * scripting facility here...
     */
    write(dev->fd, ipaddress, strlen(ipaddress));
    write(dev->fd, "\n", 1);
    return 0;
}

static int VMNET_Init(DEVBLK *dev, int argc, char *argv[])
{
U16             xdevnum;                /* Pair device devnum        */
DEVBLK          *xdev;                  /* Pair device               */
int rc;
U16 lcss;

    dev->devtype = 0x3088;

    /* parameters for network CTC are:
     *    devnum of the other CTC device of the pair
     *    ipaddress
     *    vmnet command line
     *
     * CTC adapters are used in pairs, one for READ, one for WRITE.
     * The vmnet is only initialised when both are initialised.
     */
    if (argc < 3) {
        logmsg(_("HHCCT027E %4.4X: Not enough parameters\n"), dev->devnum);
        return -1;
    }
    rc=parse_single_devnum(argv[0],&lcss,&xdevnum);
    if (rc<0)
    {
        logmsg(_("HHCCT028E %d:%4.4X: Bad device number '%s'\n"),
                  SSID_TO_LCSS(dev->ssid), dev->devnum, argv[0]);
        return -1;
    }
    xdev = find_device_by_devnum(lcss,xdevnum);
    if (xdev != NULL) {
        if (start_vmnet(dev, xdev, argc - 1, &argv[1]))
            return -1;
    }
    strcpy(dev->filename, "vmnet");

    /* Set the control unit type */
    /* Linux/390 currently only supports 3088 model 2 CTCA and ESCON */
    dev->ctctype = CTC_VMNET;

    SetSIDInfo( dev, 0x3088, 0x08, 0x3088, 0x01 );

    /* Initialize the device dependent fields */
    dev->ctcpos = 0;
    dev->ctcrem = 0;

    /* Set length of buffer */
    /* This size guarantees we can write a full iobuf of 65536
     * as a SLIP packet in a single write.  Probably overkill... */
    dev->bufsize = 65536 * 2 + 1;
    return 0;
}

static int VMNET_Write(DEVBLK *dev, BYTE *iobuf, U16 count, BYTE *unitstat)
{
int blklen = (iobuf[0]<<8) | iobuf[1];
int pktlen;
BYTE *p = iobuf + 2;
BYTE *buffer = dev->buf;
int len = 0, rem;

    if (count < blklen) {
        logmsg (_("HHCCT029E %4.4X: bad block length: %d < %d\n"),
                dev->devnum, count, blklen);
        blklen = count;
    }
    while (p < iobuf + blklen) {
        pktlen = (p[0]<<8) | p[1];

        rem = iobuf + blklen - p;

        if (rem < pktlen) {
            logmsg (_("HHCCT030E %4.4X: bad packet length: %d < %d\n"),
                    dev->devnum, rem, pktlen);
            pktlen = rem;
        }
        if (pktlen < 6) {
        logmsg (_("HHCCT031E %4.4X: bad packet length: %d < 6\n"),
                    dev->devnum, pktlen);
            pktlen = 6;
        }

        pktlen -= 6;
        p += 6;

        while (pktlen--) {
            switch (*p) {
            case SLIP_END:
                buffer[len++] = SLIP_ESC;
                buffer[len++] = SLIP_ESC_END;
                break;
            case SLIP_ESC:
                buffer[len++] = SLIP_ESC;
                buffer[len++] = SLIP_ESC_ESC;
                break;
            default:
                buffer[len++] = *p;
                break;
            }
            p++;
        }
        buffer[len++] = SLIP_END;
        write(dev->fd, buffer, len);   /* should check error conditions? */
        len = 0;
    }

    *unitstat = CSW_CE | CSW_DE;

    return count;
}

static int bufgetc(DEVBLK *dev, int blocking)
{
BYTE *bufp = dev->buf + dev->ctcpos, *bufend = bufp + dev->ctcrem;
int n;

    if (bufp >= bufend) {
        if (blocking == 0) return -1;
        do {
            n = read(dev->fd, dev->buf, dev->bufsize);
            if (n <= 0) {
                if (n == 0) {
                    /* VMnet died on us. */
                    logmsg (_("HHCCT032E %4.4X: Error: EOF on read, "
                              "CTC network down\n"),
                            dev->devnum);
                    /* -2 will cause an error status to be set */
                    return -2;
                }
                if( n == EINTR )
                    return -3;
                logmsg (_("HHCCT033E %4.4X: Error: read: %s\n"),
                        dev->devnum, strerror(errno));
                SLEEP(2);
            }
        } while (n <= 0);
        dev->ctcrem = n;
        bufend = &dev->buf[n];
        dev->ctclastpos = dev->ctclastrem = dev->ctcpos = 0;
        bufp = dev->buf;
    }

    dev->ctcpos++;
    dev->ctcrem--;

    return *bufp;
}

static void setblkheader(BYTE *iobuf, int buflen)
{
    iobuf[0] = (buflen >> 8) & 0xFF;
    iobuf[1] = buflen & 0xFF;
}

static void setpktheader(BYTE *iobuf, int packetpos, int packetlen)
{
    iobuf[packetpos] = (packetlen >> 8) & 0xFF;
    iobuf[packetpos+1] = packetlen & 0xFF;
    iobuf[packetpos+2] = 0x08;
    iobuf[packetpos+3] = 0;
    iobuf[packetpos+4] = 0;
    iobuf[packetpos+5] = 0;
}

/* read data from the CTC connection.
 * If a packet overflows the iobuf or the read buffer runs out, there are
 * 2 possibilities:
 * - block has single packet: continue reading packet, drop bytes,
 *   then return truncated packet.
 * - block has multiple packets: back up on last packet and return
 *   what we have.  Do this last packet in the next IO.
 */
static int VMNET_Read(DEVBLK *dev, BYTE *iobuf, U16 count, BYTE *unitstat)
{
int             c;                      /* next byte to process      */
int             len = 8;                /* length of block           */
int             lastlen = 2;            /* block length at last pckt */

    dev->ctclastpos = dev->ctcpos;
    dev->ctclastrem = dev->ctcrem;

    while (1) {
        c = bufgetc(dev, lastlen == 2);
        if (c < 0) {
            if(c == -3)
                return 0;
            /* End of input buffer.  Return what we have. */

            setblkheader (iobuf, lastlen);

            dev->ctcpos = dev->ctclastpos;
            dev->ctcrem = dev->ctclastrem;

            *unitstat = CSW_CE | CSW_DE | (c == -2 ? CSW_UX : 0);

            return lastlen;
        }
        switch (c) {
        case SLIP_END:
            if (len > 8) {
                /* End of packet.  Set up for next. */

                setpktheader (iobuf, lastlen, len-lastlen);

                dev->ctclastpos = dev->ctcpos;
                dev->ctclastrem = dev->ctcrem;
                lastlen = len;

                len += 6;
            }
            break;
        case SLIP_ESC:
            c = bufgetc(dev, lastlen == 2);
            if (c < 0) {
                if(c == -3)
                    return 0;
                /* End of input buffer.  Return what we have. */

                setblkheader (iobuf, lastlen);

                dev->ctcpos = dev->ctclastpos;
                dev->ctcrem = dev->ctclastrem;

                *unitstat = CSW_CE | CSW_DE | (c == -2 ? CSW_UX : 0);

                return lastlen;
            }
            switch (c) {
            case SLIP_ESC_END:
                c = SLIP_END;
                break;
            case SLIP_ESC_ESC:
                c = SLIP_ESC;
                break;
            }
            /* FALLTHRU */
        default:
            if (len < count) {
                iobuf[len++] = c;
            } else if (lastlen > 2) {
                /* IO buffer is full and we have data to return */

                setblkheader (iobuf, lastlen);

                dev->ctcpos = dev->ctclastpos;
                dev->ctcrem = dev->ctclastrem;

                *unitstat = CSW_CE | CSW_DE | (c == -2 ? CSW_UX : 0);

                return lastlen;
            } /* else truncate end of very large single packet... */
        }
    }
}
/*-------------------------------------------------------------------*/
/* End of VMNET functions written by Willem Konynenberg              */
/*-------------------------------------------------------------------*/

// ====================================================================
// Support Functions
// ====================================================================

// ---------------------------------------------------------------------
// ParseMAC
// ---------------------------------------------------------------------
//
// Parse a string containing a MAC (hardware) address and return the
// binary equivalent.
//
// Input:
//      pszMACAddr   Pointer to string containing a MAC Address in the
//                   format "xx-xx-xx-xx-xx-xx" or "xx:xx:xx:xx:xx:xx".
//
// Output:
//      pbMACAddr    Pointer to a BYTE array to receive the MAC Address
//                   that MUST be at least sizeof(MAC) bytes long.
//
// Returns:
//      0 on success, -1 otherwise
//

int  ParseMAC( char* pszMACAddr, BYTE* pbMACAddr )
{
    char    work[((sizeof(MAC)*3)-0)];
    BYTE    sep;
    int       x;
    unsigned  i;

    if (strlen(pszMACAddr) != ((sizeof(MAC)*3)-1)
        || (sizeof(MAC) > 1 &&
            *(pszMACAddr+2) != '-' &&
            *(pszMACAddr+2) != ':')
    )
    {
        errno = EINVAL;
        return -1;
    }

    strncpy(work,pszMACAddr,((sizeof(MAC)*3)-1));
    work[((sizeof(MAC)*3)-1)] = sep = *(pszMACAddr+2);

    for (i=0; i < sizeof(MAC); i++)
    {
        if
        (0
            || !isxdigit(work[(i*3)+0])
            || !isxdigit(work[(i*3)+1])
            ||  sep  !=  work[(i*3)+2]
        )
        {
            errno = EINVAL;
            return -1;
        }

        work[(i*3)+2] = 0;
        sscanf(&work[(i*3)+0],"%x",&x);
        *(pbMACAddr+i) = x;
    }

    return 0;
}

// ---------------------------------------------------------------------
// packet_trace
// ---------------------------------------------------------------------
//
// Subroutine to trace the contents of a buffer
//

void  packet_trace( BYTE* pAddr, int iLen )
{
    int           offset;
    unsigned int  i;
    unsigned char c = '\0';
    unsigned char e = '\0';
    unsigned char print_chars[17];

    for( offset = 0; offset < iLen; )
    {
        memset( print_chars, 0, sizeof( print_chars ) );

        logmsg( "+%4.4X  ", offset );

        for( i = 0; i < 16; i++ )
        {
            c = *pAddr++;

            if( offset < iLen )
            {
                logmsg("%2.2X", c);

                print_chars[i] = '.';
                e = guest_to_host( c );

                if( isprint( e ) )
                    print_chars[i] = e;
                if( isprint( c ) )
                    print_chars[i] = c;
            }
            else
            {
                logmsg( "  " );
            }

            offset++;
            if( ( offset & 3 ) == 0 )
            {
                logmsg( " " );
            }
        }

        logmsg( " %s\n", print_chars );
    }
}
