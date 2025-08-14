/* FAKETAPE.C   (c) Copyright "Fish" (David B. Trout), 2009          */
/*              Hercules Tape Device Handler for FAKETAPE            */

/*-------------------------------------------------------------------*/
/* This module contains the FAKETAPE emulated tape format support.   */
/*                                                                   */
/* The subroutines in this module are called by the general tape     */
/* device handler (tapedev.c) when the tape format is FAKETAPE.      */
/*                                                                   */
/* Messages issued by this module are prefixed HHCTA5nn              */
/*-------------------------------------------------------------------*/

/*-------------------------------------------------------------------*/
/* Reference information:                                            */
/* FSIMS100 Faketape manual                                          */
/*-------------------------------------------------------------------*/

#include "hstdinc.h"
#include "hercules.h"  /* need Hercules control blocks               */
#include "tapedev.h"   /* Main tape handler header file              */

//#define  ENABLE_TRACING_STMTS     // (Fish: DEBUGGING)

#ifdef ENABLE_TRACING_STMTS
  #if !defined(DEBUG)
    #warning DEBUG required for ENABLE_TRACING_STMTS
  #endif
  // (TRACE, ASSERT, and VERIFY macros are #defined in hmacros.h)
#else
  #undef  TRACE
  #undef  ASSERT
  #undef  VERIFY
  #define TRACE       1 ? ((void)0) : logmsg
  #define ASSERT(a)
  #define VERIFY(a)   ((void)(a))
#endif

/*********************************************************************/
/* START OF ORIGINAL AWS FUNCTIONS   (ISW Additions)                 */
/*********************************************************************/

/*-------------------------------------------------------------------*/
/* Close a FAKETAPE format file                                      */
/* New Function added by ISW for consistency with other medias       */
/*-------------------------------------------------------------------*/
void close_faketape (DEVBLK *dev)
{
    if(dev->fd>=0)
    {
        logmsg (_("HHCTA501I %4.4X: FakeTape %s closed\n"),
                dev->devnum, dev->filename);
        close(dev->fd);
    }
    strcpy(dev->filename, TAPE_UNLOADED);
    dev->fd=-1;
    dev->blockid = 0;
    dev->fenced = 0;
    return;
}

/*-------------------------------------------------------------------*/
/* Rewinds a FAKETAPE format file                                    */
/* New Function added by ISW for consistency with other medias       */
/*-------------------------------------------------------------------*/
int rewind_faketape (DEVBLK *dev,BYTE *unitstat,BYTE code)
{
    off_t rcoff;
    rcoff=lseek(dev->fd,0,SEEK_SET);
    if(rcoff<0)
    {
        build_senseX(TAPE_BSENSE_REWINDFAILED,dev,unitstat,code);
        return -1;
    }
    dev->nxtblkpos=0;
    dev->prvblkpos=-1;
    dev->curfilen=1;
    dev->blockid=0;
    dev->fenced = 0;
    return 0;
}

/*-------------------------------------------------------------------*/
/* Determines if a FAKETAPE has passed a virtual EOT marker          */
/* New Function added by ISW for consistency with other medias       */
/*-------------------------------------------------------------------*/
int passedeot_faketape (DEVBLK *dev)
{
    if(dev->nxtblkpos==0)
    {
        dev->eotwarning = 0;
        return 0;
    }
    if(dev->tdparms.maxsize==0)
    {
        dev->eotwarning = 0;
        return 0;
    }
    if(dev->nxtblkpos+dev->eotmargin > dev->tdparms.maxsize)
    {
        dev->eotwarning = 1;
        return 1;
    }
    dev->eotwarning = 0;
    return 0;
}

/*********************************************************************/
/* START OF ORIGINAL RB AWS FUNCTIONS                                */
/*********************************************************************/

/*-------------------------------------------------------------------*/
/* Open a FAKETAPE format file                                       */
/*                                                                   */
/* If successful, the file descriptor is stored in the device block  */
/* and the return value is zero.  Otherwise the return value is -1.  */
/*-------------------------------------------------------------------*/
int open_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc = -1;                /* Return code               */
char            pathname[MAX_PATH];     /* file path in host format  */

    /* Check for no tape in drive */
    if (!strcmp (dev->filename, TAPE_UNLOADED))
    {
        build_senseX(TAPE_BSENSE_TAPEUNLOADED,dev,unitstat,code);
        return -1;
    }

    /* Open the FAKETAPE file */
    hostpath(pathname, dev->filename, sizeof(pathname));
    if(!dev->tdparms.logical_readonly)
    {
        rc = hopen(pathname, O_RDWR | O_BINARY);
    }

    /* If file is read-only, attempt to open again */
    if (dev->tdparms.logical_readonly || (rc < 0 && (EROFS == errno || EACCES == errno)))
    {
        dev->readonly = 1;
        rc = hopen(pathname, O_RDONLY | O_BINARY);
    }

    /* Check for successful open */
    if (rc < 0)
    {
        logmsg (_("HHCTA502E %4.4X: Error opening %s: %s\n"),
                dev->devnum, dev->filename, strerror(errno));

        strcpy(dev->filename, TAPE_UNLOADED);
        build_senseX(TAPE_BSENSE_TAPELOADFAIL,dev,unitstat,code);
        return -1;
    }

    /* Store the file descriptor in the device block */
    dev->fd = rc;
    rc=rewind_faketape(dev,unitstat,code);
    return rc;

} /* end function open_faketape */

/*-------------------------------------------------------------------*/
/* Read a FAKETAPE block header                                      */
/*                                                                   */
/* If successful, return value is zero, and prvblkl and curblkl are  */
/* set to the previous and current block lengths respectively.       */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/* and prvblkl and curblkl are undefined. Either or both of prvblkl  */
/* and/or curblkl may be NULL.                                       */
/*-------------------------------------------------------------------*/
int readhdr_faketape (DEVBLK *dev, off_t blkpos,
                      U16* pprvblkl, U16* pcurblkl,
                      BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
off_t           rcoff;                  /* Return code from lseek()  */
FAKETAPE_BLKHDR fakehdr;                /* FakeTape block header     */
char            sblklen[5];             /* work for converting hdr   */
int             prvblkl;                /* Previous block length     */
int             curblkl;                /* Current block length      */
int             xorblkl;                /* XOR check of block lens   */

    /* Reposition file to the requested block header */
    rcoff = lseek (dev->fd, blkpos, SEEK_SET);
    if (rcoff < 0)
    {
        /* Handle seek error condition */
        logmsg (_("HHCTA503E %4.4X: Error seeking to offset "I64_FMTX" "
                "in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_LOCATEERR,dev,unitstat,code);
        return -1;
    }

    /* Read the 12-ASCII-hex-character block header */
    rc = read (dev->fd, &fakehdr, sizeof(FAKETAPE_BLKHDR));

    /* Handle read error condition */
    if (rc < 0)
    {
        logmsg (_("HHCTA504E %4.4X: Error reading block header "
                "at offset "I64_FMTX" in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_READFAIL,dev,unitstat,code);
        return -1;
    }

    /* Handle end of file (uninitialized tape) condition */
    if (rc == 0)
    {
        logmsg (_("HHCTA505E %4.4X: End of file (end of tape) "
                "at offset "I64_FMTX" in file %s\n"),
                dev->devnum, blkpos, dev->filename);

        /* Set unit exception with tape indicate (end of tape) */
        build_senseX(TAPE_BSENSE_EMPTYTAPE,dev,unitstat,code);
        return -1;
    }

    /* Handle end of file within block header */
    if (rc < (int)sizeof(FAKETAPE_BLKHDR))
    {
        logmsg (_("HHCTA506E %4.4X: Unexpected end of file in block header "
                "at offset "I64_FMTX" in file %s\n"),
                dev->devnum, blkpos, dev->filename);

        build_senseX(TAPE_BSENSE_BLOCKSHORT,dev,unitstat,code);
        return -1;
    }

    /* Convert the ASCII-hex-character block lengths to binary */
    strncpy( sblklen, fakehdr.sprvblkl, 4 ); sblklen[4] = 0; sscanf( sblklen, "%x", &prvblkl );
    strncpy( sblklen, fakehdr.scurblkl, 4 ); sblklen[4] = 0; sscanf( sblklen, "%x", &curblkl );
    strncpy( sblklen, fakehdr.sxorblkl, 4 ); sblklen[4] = 0; sscanf( sblklen, "%x", &xorblkl );

    /* Verify header integrity using the XOR header field */
    if ( (prvblkl ^ curblkl) != xorblkl )
    {
        logmsg (_("HHCTA507E %4.4X: Block header damage "
                "at offset "I64_FMTX" in file %s\n"),
                dev->devnum, blkpos, dev->filename);

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_READFAIL,dev,unitstat,code);
        return -1;
    }

    /* Return the converted value(s) to the caller */
    if (pprvblkl) *pprvblkl = prvblkl;
    if (pcurblkl) *pcurblkl = curblkl;

    /* Successful return */
    return 0;

} /* end function readhdr_faketape */

/*-------------------------------------------------------------------*/
/* Read a block from a FAKETAPE format file                          */
/*                                                                   */
/* If successful, return value is block length read.                 */
/* If a tapemark was read, the return value is zero, and the         */
/* current file number in the device block is incremented.           */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int read_faketape (DEVBLK *dev, BYTE *buf, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
off_t           blkpos;                 /* Offset of block header    */
U16             curblkl;                /* Current block length      */

    /* Initialize current block position */
    blkpos = dev->nxtblkpos;

    /* Read the block header to obtain the current block length */
    rc = readhdr_faketape (dev, blkpos, NULL, &curblkl, unitstat,code);
    if (rc < 0) return -1; /* (error message already issued) */
    ASSERT( curblkl >= 0 );

    /* Calculate the offset of the next block header */
    blkpos += sizeof(FAKETAPE_BLKHDR) + curblkl;

#if 0
    /*BHE following code will never be true!!*/

    /* Check that block length will not exceed buffer size */
    if (curblkl > MAX_BLKLEN)
    {
        logmsg (_("HHCTA508E %4.4X: Block length exceeds %d "
                "at offset "I64_FMTX" in file %s\n"),
                dev->devnum, (int)MAX_BLKLEN, blkpos, dev->filename);

        /* Set unit check with data check */
        build_senseX(TAPE_BSENSE_READFAIL,dev,unitstat,code);
        return -1;
    }
#endif

    /* If not a tapemark, read the data block */
    if (curblkl > 0)
    {
        rc = read (dev->fd, buf, curblkl);

        /* Handle read error condition */
        if (rc < 0)
        {
            logmsg (_("HHCTA510E %4.4X: Error reading data block "
                    "at offset "I64_FMTX" in file %s: %s\n"),
                    dev->devnum, blkpos, dev->filename, strerror(errno));

            /* Set unit check with equipment check */
            build_senseX(TAPE_BSENSE_READFAIL,dev,unitstat,code);
            return -1;
        }

        /* Handle end of file within data block */
        if (rc < curblkl)
        {
            logmsg (_("HHCTA511E %4.4X: Unexpected end of file in data block "
                    "at offset "I64_FMTX" in file %s\n"),
                    dev->devnum, blkpos, dev->filename);

            /* Set unit check with data check and partial record */
            build_senseX(TAPE_BSENSE_BLOCKSHORT,dev,unitstat,code);
            return -1;
        }
    }

    /* Calculate the offsets of the next and previous blocks */
    dev->prvblkpos = dev->nxtblkpos;
    dev->nxtblkpos = blkpos;

    /* Increment the block number */
    dev->blockid++;

    /* Increment file number and return zero if tapemark was read */
    if (curblkl == 0)
    {
        dev->curfilen++;
        return 0; /* UX will be set by caller */
    }

    /* Return block length */
    return curblkl;

} /* end function read_faketape */

/*-------------------------------------------------------------------*/
/* Write a FAKETAPE block header                                     */
/*                                                                   */
/* If successful, return value is zero.                              */
/* If error, return value is -1 and unitstat is set to CE+DE+UC.     */
/*-------------------------------------------------------------------*/

int  writehdr_faketape (DEVBLK *dev, off_t blkpos,
                        U16 prvblkl, U16 curblkl,
                        BYTE *unitstat, BYTE code)
{
int             rc;                     /* Return code               */
off_t           rcoff;                  /* Return code from lseek()  */
FAKETAPE_BLKHDR fakehdr;                /* FAKETAPE block header     */
char            sblklen[5];             /* work buffer               */

    /* Position file to where block header is to go */
    rcoff = lseek (dev->fd, blkpos, SEEK_SET);
    if (rcoff < 0)
    {
        /* Handle seek error condition */
        logmsg (_("HHCTA512E %4.4X: Error seeking to offset "I64_FMTX" "
                "in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_LOCATEERR,dev,unitstat,code);
        return -1;
    }

    /* Build the 12-ASCII-hex-character block header */
    snprintf( sblklen, sizeof(sblklen), "%4.4X", prvblkl );
    memcpy( fakehdr.sprvblkl, sblklen, sizeof(fakehdr.sprvblkl) );
    snprintf( sblklen, sizeof(sblklen), "%4.4X", curblkl );
    memcpy( fakehdr.scurblkl, sblklen, sizeof(fakehdr.scurblkl) );
    snprintf( sblklen, sizeof(sblklen), "%4.4X", prvblkl ^ curblkl );
    memcpy( fakehdr.sxorblkl, sblklen, sizeof(fakehdr.sxorblkl) );

    /* Write the block header */
    rc = write (dev->fd, &fakehdr, sizeof(FAKETAPE_BLKHDR));
    if (rc < (int)sizeof(FAKETAPE_BLKHDR))
    {
        if(errno==ENOSPC)
        {
            /* Disk FULL */
            build_senseX(TAPE_BSENSE_ENDOFTAPE,dev,unitstat,code);
            logmsg (_("HHCTA513E %4.4X: Media full condition reached "
                    "at offset "I64_FMTX" in file %s\n"),
                    dev->devnum, blkpos, dev->filename);
            return -1;
        }
        /* Handle write error condition */
        logmsg (_("HHCTA514E %4.4X: Error writing block header "
                "at offset "I64_FMTX" in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_WRITEFAIL,dev,unitstat,code);
        return -1;
    }

    return 0;

} /* end function writehdr_faketape */

/*-------------------------------------------------------------------*/
/* Write a block to a FAKETAPE format file                           */
/*                                                                   */
/* If successful, return value is zero.                              */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int write_faketape (DEVBLK *dev, BYTE *buf, U16 blklen,
                        BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
off_t           rcoff;                  /* Return code from lseek()  */
off_t           blkpos;                 /* Offset of block header    */
U16             prvblkl;                /* Length of previous block  */

    /* Initialize current block position and previous block length */
    blkpos = dev->nxtblkpos;
    prvblkl = 0;

    /* Determine previous block length if not at start of tape */
    if (dev->nxtblkpos > 0)
    {
        /* Retrieve the previous block length */
        rc = readhdr_faketape (dev, dev->prvblkpos, NULL, &prvblkl, unitstat,code);
        if (rc < 0) return -1;

        /* Recalculate the offset of the next block */
        blkpos = dev->prvblkpos + sizeof(FAKETAPE_BLKHDR) + prvblkl;
    }

    /* Reposition file to the new block header */
    rcoff = lseek (dev->fd, blkpos, SEEK_SET);
    if (rcoff < 0)
    {
        /* Handle seek error condition */
        logmsg (_("HHCTA515E %4.4X: Error seeking to offset "I64_FMTX" "
                "in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_LOCATEERR,dev,unitstat,code);
        return -1;
    }
    /* ISW: Determine if we are passed maxsize */
    if(dev->tdparms.maxsize>0)
    {
        if((off_t)(dev->nxtblkpos+blklen+sizeof(FAKETAPE_BLKHDR)) > dev->tdparms.maxsize)
        {
            build_senseX(TAPE_BSENSE_ENDOFTAPE,dev,unitstat,code);
            return -1;
        }
    }
    /* ISW: End of virtual physical EOT determination */

    /* Write the block header */
    rc = writehdr_faketape (dev, rcoff, prvblkl, blklen, unitstat, code);
    if (rc < 0) return -1; /* (error message already issued) */

    /* Calculate the offsets of the next and previous blocks */
    dev->nxtblkpos = blkpos + sizeof(FAKETAPE_BLKHDR) + blklen;
    dev->prvblkpos = blkpos;

    /* Write the data block */
    rc = write (dev->fd, buf, blklen);
    if (rc < blklen)
    {
        if(errno==ENOSPC)
        {
            /* Disk FULL */
            build_senseX(TAPE_BSENSE_ENDOFTAPE,dev,unitstat,code);
            logmsg (_("HHCTA516E %4.4X: Media full condition reached "
                    "at offset "I64_FMTX" in file %s\n"),
                    dev->devnum, blkpos, dev->filename);
            return -1;
        }
        /* Handle write error condition */
        logmsg (_("HHCTA517E %4.4X: Error writing data block "
                "at offset "I64_FMTX" in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_WRITEFAIL,dev,unitstat,code);
        return -1;
    }

    /* Increment the block number */
    dev->blockid++;

    /* Set new physical EOF */
    do rc = ftruncate( dev->fd, dev->nxtblkpos );
    while (EINTR == rc);

    if (rc != 0)
    {
        /* Handle write error condition */
        logmsg (_("HHCTA518E %4.4X: Error writing data block "
                "at offset "I64_FMTX" in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_WRITEFAIL,dev,unitstat,code);
        return -1;
    }

    /* Return normal status */
    return 0;

} /* end function write_faketape */

/*-------------------------------------------------------------------*/
/* Write a tapemark to a FAKETAPE format file                        */
/*                                                                   */
/* If successful, return value is zero.                              */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int write_fakemark (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
off_t           rcoff;                  /* Return code from lseek()  */
off_t           blkpos;                 /* Offset of block header    */
U16             prvblkl;                /* Length of previous block  */

    /* Initialize current block position and previous block length */
    blkpos = dev->nxtblkpos;
    prvblkl = 0;

    /* Determine previous block length if not at start of tape */
    if (dev->nxtblkpos > 0)
    {
        /* Retrieve the previous block length */
        rc = readhdr_faketape (dev, dev->prvblkpos, NULL, &prvblkl, unitstat,code);
        if (rc < 0) return -1;

        /* Recalculate the offset of the next block */
        blkpos = dev->prvblkpos + sizeof(FAKETAPE_BLKHDR) + prvblkl;
    }

    /* Reposition file to the new block header */
    rcoff = lseek (dev->fd, blkpos, SEEK_SET);
    if (rcoff < 0)
    {
        /* Handle seek error condition */
        logmsg (_("HHCTA519E %4.4X: Error seeking to offset "I64_FMTX" "
                "in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        build_senseX(TAPE_BSENSE_LOCATEERR,dev,unitstat,code);
        return -1;
    }
    /* ISW: Determine if we are passed maxsize */
    if(dev->tdparms.maxsize>0)
    {
        if((off_t)(dev->nxtblkpos+sizeof(FAKETAPE_BLKHDR)) > dev->tdparms.maxsize)
        {
            build_senseX(TAPE_BSENSE_ENDOFTAPE,dev,unitstat,code);
            return -1;
        }
    }
    /* ISW: End of virtual physical EOT determination */

    /* Write the block header */
    rc = writehdr_faketape (dev, rcoff, prvblkl, 0, unitstat, code);
    if (rc < 0) return -1; /* (error message already issued) */

    /* Increment the block number */
    dev->blockid++;

    /* Calculate the offsets of the next and previous blocks */
    dev->nxtblkpos = blkpos + sizeof(FAKETAPE_BLKHDR);
    dev->prvblkpos = blkpos;

    /* Set new physical EOF */
    do rc = ftruncate( dev->fd, dev->nxtblkpos );
    while (EINTR == rc);

    if (rc != 0)
    {
        /* Handle write error condition */
        logmsg (_("HHCTA520E %4.4X: Error writing tape mark "
                "at offset "I64_FMTX" in file %s: %s\n"),
                dev->devnum, blkpos, dev->filename, strerror(errno));

        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_WRITEFAIL,dev,unitstat,code);
        return -1;
    }

    /* Return normal status */
    return 0;

} /* end function write_fakemark */

/*-------------------------------------------------------------------*/
/* Synchronize a FAKETAPE format file  (i.e. flush buffers to disk)  */
/*                                                                   */
/* If successful, return value is zero.                              */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int sync_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
    /* Unit check if tape is write-protected */
    if (dev->readonly)
    {
        build_senseX(TAPE_BSENSE_WRITEPROTECT,dev,unitstat,code);
        return -1;
    }

    /* Perform sync. Return error on failure. */
    if (fdatasync( dev->fd ) < 0)
    {
        /* Log the error */
        logmsg (_("HHCTA521E %4.4X: Sync error on file %s: %s\n"),
            dev->devnum, dev->filename, strerror(errno));
        /* Set unit check with equipment check */
        build_senseX(TAPE_BSENSE_WRITEFAIL,dev,unitstat,code);
        return -1;
    }

    /* Return normal status */
    return 0;

} /* end function sync_faketape */

/*-------------------------------------------------------------------*/
/* Forward space over next block of a FAKETAPE format file           */
/*                                                                   */
/* If successful, return value is the length of the block skipped.   */
/* If the block skipped was a tapemark, the return value is zero,    */
/* and the current file number in the device block is incremented.   */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int fsb_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
off_t           blkpos;                 /* Offset of block header    */
U16             blklen;                 /* Block length              */

    /* Initialize current block position */
    blkpos = dev->nxtblkpos;

    /* Read the block header to obtain the current block length */
    rc = readhdr_faketape (dev, blkpos, NULL, &blklen, unitstat,code);
    if (rc < 0) return -1; /* (error message already issued) */

    /* Calculate the offset of the next block */
    blkpos += sizeof(FAKETAPE_BLKHDR) + blklen;

    /* Calculate the offsets of the next and previous blocks */
    dev->prvblkpos = dev->nxtblkpos;
    dev->nxtblkpos = blkpos;

    /* Increment current file number if tapemark was skipped */
    if (blklen == 0)
        dev->curfilen++;

    /* Increment the block number */
    dev->blockid++;

    /* Return block length or zero if tapemark */
    return blklen;

} /* end function fsb_faketape */

/*-------------------------------------------------------------------*/
/* Backspace to previous block of a FAKETAPE format file             */
/*                                                                   */
/* If successful, return value is the length of the block.           */
/* If the block is a tapemark, the return value is zero,             */
/* and the current file number in the device block is decremented.   */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int bsb_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */
U16             curblkl;                /* Length of current block   */
U16             prvblkl;                /* Length of previous block  */
off_t           blkpos;                 /* Offset of block header    */

    /* Unit check if already at start of tape */
    if (dev->nxtblkpos == 0)
    {
        build_senseX(TAPE_BSENSE_LOADPTERR,dev,unitstat,code);
        return -1;
    }

    /* Backspace to previous block position */
    blkpos = dev->prvblkpos;

    /* Read the block header to obtain the block lengths */
    rc = readhdr_faketape (dev, blkpos, &prvblkl, &curblkl, unitstat,code);
    if (rc < 0) return -1; /* (error message already issued) */

    /* Calculate the offset of the previous block */
    dev->prvblkpos = blkpos - sizeof(FAKETAPE_BLKHDR) - prvblkl;
    dev->nxtblkpos = blkpos;

    /* Decrement current file number if backspaced over tapemark */
    if (curblkl == 0)
        dev->curfilen--;

    /* Decrement the block number */
    dev->blockid--;

    /* Return block length or zero if tapemark */
    return curblkl;

} /* end function bsb_faketape */

/*-------------------------------------------------------------------*/
/* Forward space to next logical file of a FAKETAPE format file      */
/*                                                                   */
/* For FAKETAPE files, the forward space file operation is achieved  */
/* by forward spacing blocks until positioned just after a tapemark. */
/*                                                                   */
/* If successful, return value is zero, and the current file number  */
/* in the device block is incremented by fsb_faketape.               */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int fsf_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */

    while (1)
    {
        /* Forward space over next block */
        rc = fsb_faketape (dev, unitstat,code);
        if (rc < 0) return -1; /* (error message already issued) */

        /* Exit loop if spaced over a tapemark */
        if (rc == 0) break;

    } /* end while */

    /* Return normal status */
    return 0;

} /* end function fsf_faketape */

/*-------------------------------------------------------------------*/
/* Backspace to previous logical file of a FAKETAPE format file      */
/*                                                                   */
/* For FAKETAPE files, the backspace file operation is achieved      */
/* by backspacing blocks until positioned just before a tapemark     */
/* or until positioned at start of tape.                             */
/*                                                                   */
/* If successful, return value is zero, and the current file number  */
/* in the device block is decremented by bsb_faketape.               */
/* If error, return value is -1 and unitstat is set to CE+DE+UC      */
/*-------------------------------------------------------------------*/
int bsf_faketape (DEVBLK *dev, BYTE *unitstat,BYTE code)
{
int             rc;                     /* Return code               */

    while (1)
    {
        /* Exit if now at start of tape */
        if (dev->nxtblkpos == 0)
        {
            build_senseX(TAPE_BSENSE_LOADPTERR,dev,unitstat,code);
            return -1;
        }

        /* Backspace to previous block position */
        rc = bsb_faketape (dev, unitstat,code);
        if (rc < 0) return -1; /* (error message already issued) */

        /* Exit loop if backspaced over a tapemark */
        if (rc == 0) break;

    } /* end while */

    /* Return normal status */
    return 0;

} /* end function bsf_faketape */

/*********************************************************************/
/*  END OF ORIGINAL RB AWS FUNCTIONS                                 */
/*********************************************************************/
