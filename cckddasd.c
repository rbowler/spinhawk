/* CCKDDASD.C   (c) Copyright Greg Smith, 2000-2009                  */
/*              Compressed CKD Direct Access Storage Device Handler  */

/*-------------------------------------------------------------------*/
/* This module contains device functions for compressed emulated     */
/* count-key-data direct access storage devices.                     */
/*-------------------------------------------------------------------*/

#include "hstdinc.h"

#define _CCKDDASD_C_
#define _HDASD_DLL_
#include "hercules.h"
#include "devtype.h"
#include "opcode.h"

/*-------------------------------------------------------------------*/
/* Internal functions                                                */
/*-------------------------------------------------------------------*/
int     cckddasd_init(int argc, BYTE *argv[]);
int     cckddasd_term();
int     cckddasd_init_handler( DEVBLK *dev, int argc, char *argv[] );
int     cckddasd_close_device(DEVBLK *dev);
void    cckddasd_start(DEVBLK *dev);
void    cckddasd_end(DEVBLK *dev);

int     cckd_open (DEVBLK *dev, int sfx, int flags, mode_t mode);
int     cckd_close (DEVBLK *dev, int sfx);
int     cckd_read (DEVBLK *dev, int sfx, off_t off, void *buf, size_t len);
int     cckd_write (DEVBLK *dev, int sfx, off_t off, void *buf, size_t len);
int     cckd_ftruncate(DEVBLK *dev, int sfx, off_t off);
void   *cckd_malloc(DEVBLK *dev, char *id, size_t size);
void   *cckd_calloc(DEVBLK *dev, char *id, size_t n, size_t size);
void   *cckd_free(DEVBLK *dev, char *id,void *p);

int     cckd_read_track(DEVBLK *dev, int trk, BYTE *unitstat);
int     cckd_update_track(DEVBLK *dev, int trk, int off,
                         BYTE *buf, int len, BYTE *unitstat);
int     cckd_used(DEVBLK *dev);
int     cfba_read_block(DEVBLK *dev, int blkgrp, BYTE *unitstat);
int     cfba_write_block(DEVBLK *dev, int blkgrp, int off,
                         BYTE *buf, int wrlen, BYTE *unitstat);
int     cfba_used(DEVBLK *dev);
int     cckd_read_trk(DEVBLK *dev, int trk, int ra, BYTE *unitstat);
void    cckd_readahead(DEVBLK *dev, int trk);
int     cckd_readahead_scan(int *answer, int ix, int i, void *data);
void    *cckd_ra();
void    cckd_flush_cache(DEVBLK *dev);
int     cckd_flush_cache_scan(int *answer, int ix, int i, void *data);
void    cckd_flush_cache_all();
void    cckd_purge_cache(DEVBLK *dev);
int     cckd_purge_cache_scan(int *answer, int ix, int i, void *data);
void    *cckd_writer(void *arg);
int     cckd_writer_scan(int *o, int ix, int i, void *data);
off_t   cckd_get_space(DEVBLK *dev, int *size, int flags);
void    cckd_rel_space(DEVBLK *dev, off_t pos, int len, int size);
void    cckd_flush_space(DEVBLK *dev);
int     cckd_read_chdr(DEVBLK *dev);
int     cckd_write_chdr(DEVBLK *dev);
int     cckd_read_l1(DEVBLK *dev);
int     cckd_write_l1(DEVBLK *dev);
int     cckd_write_l1ent(DEVBLK *dev, int l1x);
int     cckd_read_init(DEVBLK *dev);
int     cckd_read_fsp(DEVBLK *dev);
int     cckd_write_fsp(DEVBLK *dev);
int     cckd_read_l2(DEVBLK *dev, int sfx, int l1x);
void    cckd_purge_l2(DEVBLK *dev);
int     cckd_purge_l2_scan(int *answer, int ix, int i, void *data);
int     cckd_steal_l2();
int     cckd_steal_l2_scan(int *answer, int ix, int i, void *data);
int     cckd_write_l2(DEVBLK *dev);
int     cckd_read_l2ent(DEVBLK *dev, CCKD_L2ENT *l2, int trk);
int     cckd_write_l2ent(DEVBLK *dev,   CCKD_L2ENT *l2, int trk);
int     cckd_read_trkimg(DEVBLK *dev, BYTE *buf, int trk, BYTE *unitstat);
int     cckd_write_trkimg(DEVBLK *dev, BYTE *buf, int len, int trk, int flags);
int     cckd_harden(DEVBLK *dev);
int     cckd_trklen(DEVBLK *dev, BYTE *buf);
int     cckd_null_trk(DEVBLK *dev, BYTE *buf, int trk, int nullfmt);
int     cckd_check_null_trk (DEVBLK *dev, BYTE *buf, int trk, int len);
int     cckd_cchh(DEVBLK *dev, BYTE *buf, int trk);
int     cckd_validate(DEVBLK *dev, BYTE *buf, int trk, int len);
char   *cckd_sf_name(DEVBLK *dev, int sfx);
int     cckd_sf_init(DEVBLK *dev);
int     cckd_sf_new(DEVBLK *dev);
DLL_EXPORT void   *cckd_sf_add(void *data);
DLL_EXPORT void   *cckd_sf_remove(void *data);
DLL_EXPORT void   *cckd_sf_comp(void *data);
DLL_EXPORT void   *cckd_sf_chk(void *data);
DLL_EXPORT void   *cckd_sf_stats(void *data);
int     cckd_disable_syncio(DEVBLK *dev);
void    cckd_lock_devchain(int flag);
void    cckd_unlock_devchain();
void    *cckd_gcol();
int     cckd_gc_percolate(DEVBLK *dev, unsigned int size);
int     cckd_gc_l2(DEVBLK *dev, BYTE *buf);
DEVBLK *cckd_find_device_by_devnum (U16 devnum);
BYTE   *cckd_uncompress(DEVBLK *dev, BYTE *from, int len, int maxlen, int trk);
int     cckd_uncompress_zlib(DEVBLK *dev, BYTE *to, BYTE *from, int len, int maxlen);
int     cckd_uncompress_bzip2(DEVBLK *dev, BYTE *to, BYTE *from, int len, int maxlen);
int     cckd_compress(DEVBLK *dev, BYTE **to, BYTE *from, int len, int comp, int parm);
int     cckd_compress_none(DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm);
int     cckd_compress_zlib(DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm);
int     cckd_compress_bzip2(DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm);
void    cckd_command_help();
void    cckd_command_opts();
void    cckd_command_stats();
void    cckd_command_debug();
int     cckd_command(char *op, int cmd);
DLL_EXPORT void    cckd_print_itrace();
void    cckd_trace(DEVBLK *dev, char *msg, ...);

/*-------------------------------------------------------------------*/
/* Definitions for sense data format codes and message codes         */
/*-------------------------------------------------------------------*/
#define FORMAT_0                0       /* Program or System Checks  */
#define FORMAT_1                1       /* Device Equipment Checks   */
#define FORMAT_2                2       /* 3990 Equipment Checks     */
#define FORMAT_3                3       /* 3990 Control Checks       */
#define FORMAT_4                4       /* Data Checks               */
#define FORMAT_5                5       /* Data Check + Displacement */
#define FORMAT_6                6       /* Usage Stats/Overrun Errors*/
#define FORMAT_7                7       /* Device Control Checks     */
#define FORMAT_8                8       /* Device Equipment Checks   */
#define FORMAT_9                9       /* Device Rd/Wrt/Seek Checks */
#define FORMAT_F                15      /* Cache Storage Checks      */
#define MESSAGE_0               0       /* Message 0                 */
#define MESSAGE_1               1       /* Message 1                 */
#define MESSAGE_2               2       /* Message 2                 */
#define MESSAGE_3               3       /* Message 3                 */
#define MESSAGE_4               4       /* Message 4                 */
#define MESSAGE_5               5       /* Message 5                 */
#define MESSAGE_6               6       /* Message 6                 */
#define MESSAGE_7               7       /* Message 7                 */
#define MESSAGE_8               8       /* Message 8                 */
#define MESSAGE_9               9       /* Message 9                 */
#define MESSAGE_A               10      /* Message A                 */
#define MESSAGE_B               11      /* Message B                 */
#define MESSAGE_C               12      /* Message C                 */
#define MESSAGE_D               13      /* Message D                 */
#define MESSAGE_E               14      /* Message E                 */
#define MESSAGE_F               15      /* Message F                 */

/*-------------------------------------------------------------------*/
/* Data areas                                                        */
/*-------------------------------------------------------------------*/
static  CCKD_L2ENT empty_l2[CKDDASD_NULLTRK_FMTMAX+1][256];
static  BYTE eighthexFF[] = {0xff,0xff,0xff,0xff,0xff,0xff,0xff,0xff};
DEVHND  cckddasd_device_hndinfo;
DEVHND  cfbadasd_device_hndinfo;
DLL_EXPORT CCKDBLK cckdblk;                        /* cckd global area          */

/*-------------------------------------------------------------------*/
/* CCKD global initialization                                        */
/*-------------------------------------------------------------------*/
int cckddasd_init (int argc, BYTE *argv[])
{
int             i, j;                   /* Loop indexes              */

    UNREFERENCED(argc);
    UNREFERENCED(argv);

    if (memcmp (&cckdblk.id, "CCKDBLK ", sizeof(cckdblk.id)) == 0)
        return 0;

    /* Clear the cckdblk */
    memset(&cckdblk, 0, sizeof(CCKDBLK));

    /* Initialize locks and conditions */
    memcpy (&cckdblk.id, "CCKDBLK ", sizeof(cckdblk.id));
    initialize_lock (&cckdblk.gclock);
    initialize_lock (&cckdblk.ralock);
    initialize_lock (&cckdblk.wrlock);
    initialize_lock (&cckdblk.devlock);
    initialize_condition (&cckdblk.gccond);
    initialize_condition (&cckdblk.racond);
    initialize_condition (&cckdblk.wrcond);
    initialize_condition (&cckdblk.devcond);
    initialize_condition (&cckdblk.termcond);

    /* Initialize some variables */
    cckdblk.wrprio     = 16;
    cckdblk.ranbr      = CCKD_DEFAULT_RA_SIZE;
    cckdblk.ramax      = CCKD_DEFAULT_RA;
    cckdblk.wrmax      = CCKD_DEFAULT_WRITER;
    cckdblk.gcmax      = CCKD_DEFAULT_GCOL;
    cckdblk.gcwait     = CCKD_DEFAULT_GCOLWAIT;
    cckdblk.gcparm     = CCKD_DEFAULT_GCOLPARM;
    cckdblk.readaheads = CCKD_DEFAULT_READAHEADS;
    cckdblk.freepend   = CCKD_DEFAULT_FREEPEND;
#ifdef HAVE_LIBZ
    cckdblk.comps     |= CCKD_COMPRESS_ZLIB;
#endif
#ifdef CCKD_BZIP2
    cckdblk.comps     |= CCKD_COMPRESS_BZIP2;
#endif
    cckdblk.comp       = 0xff;
    cckdblk.compparm   = -1;

    /* Initialize the readahead queue */
    cckdblk.ra1st = cckdblk.ralast = -1;
    cckdblk.rafree = 0;
    for (i = 0; i < cckdblk.ranbr; i++) cckdblk.ra[i].next = i + 1;
    cckdblk.ra[cckdblk.ranbr - 1].next = -1;

    /* Clear the empty L2 tables */
    for (i = 0; i <= CKDDASD_NULLTRK_FMTMAX; i++)
        for (j = 0; j < 256; j++)
        {
            empty_l2[i][j].pos = 0;
            empty_l2[i][j].len = empty_l2[i][j].size = i;
        }

    return 0;

} /* end function cckddasd_init */

/*-------------------------------------------------------------------*/
/* CCKD dasd global termination                                      */
/*-------------------------------------------------------------------*/
int cckddasd_term ()
{
    /* Terminate the readahead threads */
    obtain_lock (&cckdblk.ralock);
    cckdblk.ramax = 0;
    if (cckdblk.ras)
    {
        broadcast_condition (&cckdblk.racond);
        wait_condition (&cckdblk.termcond, &cckdblk.ralock);
    }
    release_lock (&cckdblk.ralock);

    /* Terminate the garbage collection threads */
    obtain_lock (&cckdblk.gclock);
    cckdblk.gcmax = 0;
    if (cckdblk.gcs)
    {
        broadcast_condition (&cckdblk.gccond);
        wait_condition (&cckdblk.termcond, &cckdblk.gclock);
    }
    release_lock (&cckdblk.gclock);

    /* Terminate the writer threads */
    obtain_lock (&cckdblk.wrlock);
    cckdblk.wrmax = 0;
    if (cckdblk.wrs)
    {
        broadcast_condition (&cckdblk.wrcond);
        wait_condition (&cckdblk.termcond, &cckdblk.wrlock);
    }
    release_lock (&cckdblk.wrlock);

    memset(&cckdblk, 0, sizeof(CCKDBLK));

    return 0;

} /* end function cckddasd_term */

/*-------------------------------------------------------------------*/
/* CKD dasd initialization                                           */
/*-------------------------------------------------------------------*/
int cckddasd_init_handler ( DEVBLK *dev, int argc, char *argv[] )
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
DEVBLK         *dev2;                   /* -> device in cckd queue   */
int             i;                      /* Counter                   */
int             fdflags;                /* File flags                */

    UNREFERENCED(argc);
    UNREFERENCED(argv);

    /* Initialize the global cckd block if necessary */
    if (memcmp (&cckdblk.id, "CCKDBLK ", sizeof(cckdblk.id)))
        cckddasd_init (0, NULL);

    /* Obtain area for cckd extension */
    dev->cckd_ext = cckd = cckd_calloc (dev, "ext", 1, sizeof(CCKDDASD_EXT));
    if (cckd == NULL)
        return -1;

    /* Initialize locks and conditions */
    initialize_lock (&cckd->iolock);
    initialize_lock (&cckd->filelock);
    initialize_condition (&cckd->iocond);

    /* Initialize some variables */
    obtain_lock (&cckd->filelock);
    cckd->l1x = cckd->sfx = cckd->l2active = -1;
    dev->cache = cckd->free1st = -1;
    cckd->fd[0] = dev->fd;
    fdflags = get_file_accmode_flags( dev->fd );
    cckd->open[0] = (fdflags & O_RDWR) ? CCKD_OPEN_RW : CCKD_OPEN_RO;
    for (i = 1; i <= CCKD_MAX_SF; i++)
    {
        cckd->fd[i] = -1;
        cckd->open[i] = CCKD_OPEN_NONE;
    }
    cckd->maxsize = sizeof(off_t) > 4 ? 0xffffffffll : 0x7fffffffll;

    /* call the chkdsk function */
    if (cckd_chkdsk (dev, 0) < 0)
        return -1;

    /* Perform initial read */
    if (cckd_read_init (dev) < 0)
        return -1;
    if (cckd->fbadasd) dev->ckdtrksz = CFBA_BLOCK_SIZE;

    /* open the shadow files */
    if (cckd_sf_init (dev) < 0)
    {
        logmsg (_("HHCCD101E %4.4X error initializing shadow files\n"), dev->devnum);
        return -1;
    }

    /* Update the device handler routines */
    if (cckd->ckddasd)
        dev->hnd = &cckddasd_device_hndinfo;
    else
        dev->hnd = &cfbadasd_device_hndinfo;
    release_lock (&cckd->filelock);

    /* Insert the device into the cckd device queue */
    cckd_lock_devchain(1);
    for (cckd = NULL, dev2 = cckdblk.dev1st; dev2; dev2 = cckd->devnext)
        cckd = dev2->cckd_ext;
    if (cckd) cckd->devnext = dev;
    else cckdblk.dev1st = dev;
    cckd_unlock_devchain();

    cckdblk.batch = dev->batch;
    if (cckdblk.batch)
    {
        cckdblk.nostress = 1;
        cckdblk.freepend = 0;
        cckdblk.linuxnull = 1;
    }

    return 0;
} /* end function cckddasd_init_handler */

/*-------------------------------------------------------------------*/
/* Close a Compressed CKD Device                                     */
/*-------------------------------------------------------------------*/
int cckddasd_close_device (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             i;                      /* Index                     */

    cckd = dev->cckd_ext;

    /* Wait for readaheads to finish */
    obtain_lock(&cckdblk.ralock);
    cckd->stopping = 1;
    while (cckd->ras)
    {
        release_lock(&cckdblk.ralock);
        usleep(1);
        obtain_lock(&cckdblk.ralock);
    }
    release_lock(&cckdblk.ralock);

    /* Flush the cache and wait for the writes to complete */
    obtain_lock (&cckd->iolock);
    cckd->stopping = 1;
    cckd_flush_cache (dev);
    while (cckd->wrpending || cckd->ioactive)
    {
        cckd->iowaiters++;
        wait_condition (&cckd->iocond, &cckd->iolock);
        cckd->iowaiters--;
        cckd_flush_cache (dev);
    }
    broadcast_condition (&cckd->iocond);
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    if (cckd->newbuf) cckd_free (dev, "newbuf", cckd->newbuf);
    release_lock (&cckd->iolock);

    /* Remove the device from the cckd queue */
    cckd_lock_devchain(1);
    if (dev == cckdblk.dev1st) cckdblk.dev1st = cckd->devnext;
    else
    {
        DEVBLK *dev2 = cckdblk.dev1st;
        CCKDDASD_EXT *cckd2 = dev2->cckd_ext;
        while (cckd2->devnext != dev)
        {
           dev2 = cckd2->devnext; cckd2 = dev2->cckd_ext;
        }
        cckd2->devnext = cckd->devnext;
    }
    cckd_unlock_devchain();

    /* harden the file */
    obtain_lock (&cckd->filelock);
    cckd_harden (dev);

    /* close the shadow files */
    for (i = 1; i <= cckd->sfn; i++)
    {
        cckd_close (dev, i);
        cckd->open[i] = 0;
    }

    /* free the level 1 tables */
    for (i = 0; i <= cckd->sfn; i++)
        cckd->l1[i] = cckd_free (dev, "l1", cckd->l1[i]);

    /* reset the device handler */
    if (cckd->ckddasd)
        dev->hnd = &ckddasd_device_hndinfo;
    else
        dev->hnd = &fbadasd_device_hndinfo;

    /* write some statistics */
    if (!dev->batch)
        cckd_sf_stats (dev);
    release_lock (&cckd->filelock);

    /* free the cckd extension */
    dev->cckd_ext= cckd_free (dev, "ext", cckd);

    if (dev->dasdsfn) free (dev->dasdsfn);
    dev->dasdsfn = NULL;

    close (dev->fd);
    dev->fd = -1;

    /* If no more devices then perform global termination */
    if (cckdblk.dev1st == NULL) cckddasd_term ();

    return 0;
} /* end function cckddasd_close_device */

/*-------------------------------------------------------------------*/
/* Compressed ckd start/resume channel program                       */
/*-------------------------------------------------------------------*/
void cckddasd_start (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
U16             devnum = 0;             /* Last active device number */
int             trk = 0;                /* Last active track         */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "start i/o file[%d] bufcur %d cache[%d]\n",
                cckd->sfn, dev->bufcur, dev->cache);

    /* Reset buffer offsets */
    dev->bufoff = 0;
    dev->bufoffhi = cckd->ckddasd ? dev->ckdtrksz : CFBA_BLOCK_SIZE;

    /* Check for merge - synchronous i/o should be disabled */
    obtain_lock(&cckd->iolock);
    if (cckd->merging)
    {
        cckd_trace (dev, "start i/o waiting for merge%s\n","");
        while (cckd->merging)
        {
            cckd->iowaiters++;
            wait_condition (&cckd->iocond, &cckd->iolock);
            cckd->iowaiters--;
        }
        dev->bufcur = dev->cache = -1;
    }
    cckd->ioactive = 1;

    cache_lock(CACHE_DEVBUF);

    if (dev->cache >= 0)
        CCKD_CACHE_GETKEY(dev->cache, devnum, trk);

    /* Check if previous active entry is still valid and not busy */
    if (dev->cache >= 0 && dev->devnum == devnum && dev->bufcur == trk
     && !(cache_getflag(CACHE_DEVBUF, dev->cache) & CCKD_CACHE_IOBUSY))
    {
        /* Make the entry active again */
        cache_setflag (CACHE_DEVBUF, dev->cache, ~0, CCKD_CACHE_ACTIVE);

        /* If the entry is pending write then change it to `updated' */
        if (cache_getflag(CACHE_DEVBUF, dev->cache) & CCKD_CACHE_WRITE)
        {
            cache_setflag (CACHE_DEVBUF, dev->cache, ~CCKD_CACHE_WRITE, CCKD_CACHE_UPDATED);
            cckd->wrpending--;
            if (cckd->iowaiters && !cckd->wrpending)
                broadcast_condition (&cckd->iocond);
        }
    }
    else
        dev->bufcur = dev->cache = -1;

    cache_unlock (CACHE_DEVBUF);

    release_lock (&cckd->iolock);

    return;

} /* end function cckddasd_start */

/*-------------------------------------------------------------------*/
/* Compressed ckd end/suspend channel program                        */
/*-------------------------------------------------------------------*/
void cckddasd_end (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */

    cckd = dev->cckd_ext;

    /* Update length if previous image was updated */
    if (dev->bufupd && dev->bufcur >= 0 && dev->cache >= 0)
    {
        dev->buflen = cckd_trklen (dev, dev->buf);
        cache_setval (CACHE_DEVBUF, dev->cache, dev->buflen);
    }

    dev->bufupd = 0;

    cckd_trace (dev, "end i/o bufcur %d cache[%d] waiters %d\n",
                dev->bufcur, dev->cache, cckd->iowaiters);

    obtain_lock (&cckd->iolock);

    cckd->ioactive = 0;

    /* Make the current entry inactive */
    if (dev->cache >= 0)
    {
        cache_lock (CACHE_DEVBUF);
        cache_setflag (CACHE_DEVBUF, dev->cache, ~CCKD_CACHE_ACTIVE, 0);
        cache_unlock (CACHE_DEVBUF);
    }

    /* Cause writers to start after first update */
    if (cckd->updated && (cckdblk.wrs == 0 || cckd->iowaiters != 0))
        cckd_flush_cache (dev);
    else if (cckd->iowaiters)
        broadcast_condition (&cckd->iocond);

    release_lock (&cckd->iolock);

} /* end function cckddasd_end */

/*-------------------------------------------------------------------*/
/* Open a cckd file                                                  */
/*                                                                   */
/* If O_CREAT is not set and mode is non-zero then the error message */
/* will be supressed.                                                */
/*-------------------------------------------------------------------*/
int cckd_open (DEVBLK *dev, int sfx, int flags, mode_t mode)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             err;                    /* 1 = issue error message   */
char            pathname[MAX_PATH];     /* file path in host format  */

    cckd = dev->cckd_ext;

    err = !((flags & O_CREAT) == 0 && mode != 0);

    if (cckd->fd[sfx] >= 0)
        cckd_close (dev, sfx);

    hostpath(pathname, cckd_sf_name (dev, sfx), sizeof(pathname));
    cckd->fd[sfx] = hopen(pathname, flags, mode);
    if (sfx == 0) dev->fd = cckd->fd[sfx];

    if (cckd->fd[sfx] >= 0)
        cckd->open[sfx] = flags & O_RDWR ? CCKD_OPEN_RW :
                          cckd->open[sfx] == CCKD_OPEN_RW ?
                          CCKD_OPEN_RD : CCKD_OPEN_RO;
    else
    {
        if (err)
        {
            logmsg (_("HHCCD130E %4.4X file[%d] %s open error: %s\n"),
                    dev->devnum, sfx, cckd_sf_name (dev, sfx),
                    strerror(errno));
            cckd_trace (dev, "file[%d] fd[%d] open %s error flags %8.8x mode %8.8x\n",
                        sfx, cckd->fd[sfx], cckd_sf_name (dev, sfx), flags, mode);
            cckd_print_itrace ();
        }
        cckd->open[sfx] = CCKD_OPEN_NONE;
    }

    cckd_trace (dev, "file[%d] fd[%d] open %s, flags %8.8x mode %8.8x\n",
                sfx, cckd->fd[sfx], cckd_sf_name (dev, sfx), flags, mode);

    return cckd->fd[sfx];

} /* end function cckd_open */

/*-------------------------------------------------------------------*/
/* Close a cckd file                                                 */
/*-------------------------------------------------------------------*/
int cckd_close (DEVBLK *dev, int sfx)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc = 0;                 /* Return code               */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "file[%d] fd[%d] close %s\n",
                sfx, cckd->fd[sfx], cckd_sf_name(dev, sfx));

    if (cckd->fd[sfx] >= 0)
        rc = close (cckd->fd[sfx]);

    if (rc < 0)
    {
        logmsg (_("HHCCD130E %4.4X file[%d] close error: %s\n"),
                dev->devnum, sfx, strerror(errno));
        cckd_print_itrace ();
    }

    cckd->fd[sfx] = -1;
    if (sfx == 0) dev->fd = -1;

    return rc;

} /* end function cckd_close */

/*-------------------------------------------------------------------*/
/* Read from a cckd file                                             */
/*-------------------------------------------------------------------*/
int cckd_read (DEVBLK *dev, int sfx, off_t off, void *buf, size_t len)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "file[%d] fd[%d] read, off 0x%" I64_FMT "x len %ld\n",
                sfx, cckd->fd[sfx], (long long)off, (long)len);

    /* Seek to specified offset */
    if (lseek (cckd->fd[sfx], off, SEEK_SET) < 0)
    {
        logmsg (_("HHCCD130E %4.4X file[%d] lseek error, offset 0x%" I64_FMT "x: %s\n"),
                dev->devnum, sfx, (long long)off, strerror(errno));
        cckd_print_itrace ();
        return -1;
    }

    /* Read the data */
    rc = read (cckd->fd[sfx], buf, len);
    if (rc < (int)len)
    {
        if (rc < 0)
            logmsg (_("HHCCD130E %4.4X file[%d] read error, offset 0x%" I64_FMT "x: %s\n"),
                    dev->devnum, sfx, (long long)off, strerror(errno));
        else
            logmsg (_("HHCCD130E %4.4X file[%d] read incomplete, offset 0x%" I64_FMT "x: "
                      "read %d expected %d\n"),
                    dev->devnum, sfx, (long long)off, rc, len);
        cckd_print_itrace ();
        return -1;
    }

    return rc;

} /* end function cckd_read */

/*-------------------------------------------------------------------*/
/* Write to a cckd file                                              */
/*-------------------------------------------------------------------*/
int cckd_write (DEVBLK *dev, int sfx, off_t off, void *buf, size_t len)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc = 0;                 /* Return code               */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "file[%d] fd[%d] write, off 0x%" I64_FMT "x len %ld\n",
                sfx, cckd->fd[sfx], (long long)off, (long)len);

    /* Seek to specified offset */
    if (lseek (cckd->fd[sfx], off, SEEK_SET) < 0)
    {
        logmsg (_("HHCCD130E %4.4X file[%d] lseek error, offset 0x%" I64_FMT "x: %s\n"),
                dev->devnum, sfx, (long long)off, strerror(errno));
        return -1;
    }

    /* Write the data */
    rc = write (cckd->fd[sfx], buf, len);
    if (rc < (int)len)
    {
        if (rc < 0)
            logmsg (_("HHCCD130E %4.4X file[%d] write error, offset 0x%" I64_FMT "x: %s\n"),
                    dev->devnum, sfx, (long long)off, strerror(errno));
        else
            logmsg (_("HHCCD130E %4.4X file[%d] write incomplete, offset 0x%" I64_FMT "x: "
                      "wrote %d expected %d\n"),
                    dev->devnum, sfx, (long long)off, rc, len);
        cckd_print_itrace ();
        return -1;
    }

    return rc;

} /* end function cckd_write */

/*-------------------------------------------------------------------*/
/* Truncate a cckd file                                              */
/*-------------------------------------------------------------------*/
int cckd_ftruncate(DEVBLK *dev, int sfx, off_t off)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "file[%d] fd[%d] ftruncate, off 0x%" I64_FMT "x\n",
                sfx, cckd->fd[sfx], (long long)off);

    /* Truncate the file */
    if (ftruncate (cckd->fd[sfx], off) < 0)
    {
        logmsg (_("HHCCD130E %4.4X file[%d] ftruncate error, offset 0x%" I64_FMT "x: %s\n"),
                dev->devnum, sfx, (long long)off, strerror(errno));
        cckd_print_itrace ();
        return -1;
    }

    return 0;

} /* end function cckd_ftruncate */

/*-------------------------------------------------------------------*/
/* malloc                                                            */
/*-------------------------------------------------------------------*/
void *cckd_malloc (DEVBLK *dev, char *id, size_t size)
{
void           *p;                      /* Pointer                   */

    p = malloc (size);
    cckd_trace (dev, "%s malloc %p len %ld\n", id, p, (long)size);

    if (p == NULL)
    {
        logmsg (_("HHCCD130E %4.4X malloc error, size %d: %s\n"),
                dev ? dev->devnum : 0, size, strerror(errno));
        cckd_print_itrace ();
    }

    return p;

} /* end function cckd_malloc */

/*-------------------------------------------------------------------*/
/* calloc                                                            */
/*-------------------------------------------------------------------*/
void *cckd_calloc (DEVBLK *dev, char *id, size_t n, size_t size)
{
void           *p;                      /* Pointer                   */

    p = calloc (n, size);
    cckd_trace (dev, "%s calloc %p len %ld\n", id, p, n*(long)size);

    if (p == NULL)
    {
        logmsg (_("HHCCD130E %4.4X calloc error, size %d: %s\n"),
                dev ? dev->devnum : 0, n*size, strerror(errno));
        cckd_print_itrace ();
    }

    return p;

} /* end function cckd_calloc */

/*-------------------------------------------------------------------*/
/* free                                                              */
/*-------------------------------------------------------------------*/
void *cckd_free (DEVBLK *dev, char *id, void *p)
{
    cckd_trace (dev, "%s free %p\n", id, p);
    if (p) free (p);
    return NULL;
} /* end function cckd_free */

/*-------------------------------------------------------------------*/
/* Compressed ckd read track image                                   */
/*-------------------------------------------------------------------*/
int cckd_read_track (DEVBLK *dev, int trk, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             len;                    /* Compressed length         */
BYTE           *newbuf;                 /* Uncompressed buffer       */
int             cache;                  /* New active cache entry    */
int             syncio;                 /* Syncio indicator          */

    cckd = dev->cckd_ext;

    /* Update length if previous image was updated */
    if (dev->bufupd && dev->bufcur >= 0 && dev->cache >= 0)
    {
        dev->buflen = cckd_trklen (dev, dev->buf);
        cache_setval (CACHE_DEVBUF, dev->cache, dev->buflen);
    }

    /* Turn off the synchronous I/O bit if trk overflow or trk 0 */
    syncio = dev->syncio_active;
    if (dev->ckdtrkof || trk == 0)
        dev->syncio_active = 0;

    /* Reset buffer offsets */
    dev->bufoff = 0;
    dev->bufoffhi = dev->ckdtrksz;

    /* Check if reading the same track image */
    if (trk == dev->bufcur && dev->cache >= 0)
    {
        /* Track image may be compressed */
        if ((dev->buf[0] & CCKD_COMPRESS_MASK) != 0
         && (dev->buf[0] & dev->comps) == 0)
        {
#if 0
            /* Return if synchronous i/o */
            if (dev->syncio_active)
            {
                cckd_trace (dev, "read  trk   %d syncio compressed\n", trk);
                cckdblk.stats_synciomisses++;
                dev->syncio_retry = 1;
                return -1;
            }
#endif
            len = cache_getval(CACHE_DEVBUF, dev->cache);
            newbuf = cckd_uncompress (dev, dev->buf, len, dev->ckdtrksz, trk);
            if (newbuf == NULL) {
                ckd_build_sense (dev, SENSE_EC, 0, 0, FORMAT_1, MESSAGE_0);
                *unitstat = CSW_CE | CSW_DE | CSW_UC;
                dev->bufcur = dev->cache = -1;
                dev->syncio_active = syncio;
                return -1;
            }
            cache_setbuf (CACHE_DEVBUF, dev->cache, newbuf, dev->ckdtrksz);
            dev->buf     = newbuf;
            dev->buflen  = cckd_trklen (dev, newbuf);
            cache_setval (CACHE_DEVBUF, dev->cache, dev->buflen);
            dev->bufsize = cache_getlen (CACHE_DEVBUF, dev->cache);
            dev->bufupd  = 0;
            cckd_trace (dev, "read  trk   %d uncompressed len %d\n",
                        trk, dev->buflen);
        }

        dev->comp = dev->buf[0] & CCKD_COMPRESS_MASK;
        if (dev->comp != 0) dev->compoff = CKDDASD_TRKHDR_SIZE;

        return 0;
    }

    cckd_trace (dev, "read  trk   %d (%s)\n", trk,
                dev->syncio_active ? "synchronous" : "asynchronous");

    /* read the new track */
    dev->bufupd = 0;
    *unitstat = 0;
    cache = cckd_read_trk (dev, trk, 0, unitstat);
    if (cache < 0)
    {
        dev->bufcur = dev->cache = -1;
        return -1;
    }

    dev->cache    = cache;
    dev->buf      = cache_getbuf (CACHE_DEVBUF, dev->cache, 0);
    dev->bufcur   = trk;
    dev->bufoff   = 0;
    dev->bufoffhi = dev->ckdtrksz;
    dev->buflen   = cache_getval (CACHE_DEVBUF, dev->cache);
    dev->bufsize  = cache_getlen (CACHE_DEVBUF, dev->cache);

    dev->comp = dev->buf[0] & CCKD_COMPRESS_MASK;
    if (dev->comp != 0) dev->compoff = CKDDASD_TRKHDR_SIZE;

    /* If the image is compressed then call ourself recursively
       to cause the image to get uncompressed */
    if (dev->comp != 0 && (dev->comp & dev->comps) == 0)
        rc = cckd_read_track (dev, trk, unitstat);
    else
        rc = 0;

    dev->syncio_active = syncio;
    return rc;
} /* end function cckd_read_track */

/*-------------------------------------------------------------------*/
/* Compressed ckd update track image                                 */
/*-------------------------------------------------------------------*/
int cckd_update_track (DEVBLK *dev, int trk, int off,
                       BYTE *buf, int len, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */

    cckd = dev->cckd_ext;

    /* Error if opened read-only */
    if (dev->ckdrdonly && cckd->sfn == 0)
    {
        ckd_build_sense (dev, SENSE_EC, SENSE1_WRI, 0,FORMAT_1, MESSAGE_0);
        *unitstat = CSW_CE | CSW_DE | CSW_UC;
        dev->bufcur = dev->cache = -1;
        return -1;
    }

    /* If the track is not current or compressed then read it.
       `dev->comps' is set to zero forcing the read routine to
       uncompress the image.                                     */
    if (trk != dev->bufcur || (dev->buf[0] & CCKD_COMPRESS_MASK) != 0)
    {
        dev->comps = 0;
        rc = (dev->hnd->read) (dev, trk, unitstat);
        if (rc < 0)
        {
            dev->bufcur = dev->cache = -1;
            return -1;
        }
    }

    /* Invalid track format if going past buffer end */
    if (off + len > dev->ckdtrksz)
    {
        ckd_build_sense (dev, 0, SENSE1_ITF, 0, 0, 0);
        *unitstat = CSW_CE | CSW_DE | CSW_UC;
        dev->bufcur = dev->cache = -1;
        return -1;
    }

    /* Copy the data into the buffer */
    if (buf && len > 0) memcpy (dev->buf + off, buf, len);

    cckd_trace (dev, "updt  trk   %d offset %d length %d\n",
                trk, off, len);

    /* Update the cache entry */
    cache_setflag (CACHE_DEVBUF, dev->cache, ~0, CCKD_CACHE_UPDATED | CCKD_CACHE_USED);
    cckd->updated = 1;

    /* Notify the shared server of the update */
    if (!dev->bufupd)
    {
        dev->bufupd = 1;
        shared_update_notify (dev, trk);
    }

    return len;

} /* end function cckd_update_track */

/*-------------------------------------------------------------------*/
/* Return used cylinders                                             */
/*-------------------------------------------------------------------*/
int cckd_used (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             l1x, l2x;               /* Lookup table indexes      */
int             sfx;                    /* Shadow file suffix        */
CCKD_L2ENT      l2;                     /* Copied level 2 entry      */

    cckd = dev->cckd_ext;
    obtain_lock (&cckd->filelock);

    /* Find the last used level 1 table entry */
    for (l1x = cckd->cdevhdr[0].numl1tab - 1; l1x > 0; l1x--)
    {
        sfx = cckd->sfn;
        while (cckd->l1[sfx][l1x] == 0xffffffff && sfx > 0) sfx--;
        if (cckd->l1[sfx][l1x]) break;
    }

    /* Find the last used level 2 table entry */
    for (l2x = 255; l2x >= 0; l2x--)
    {
        rc = cckd_read_l2ent (dev, &l2, l1x * 256 + l2x);
        if (rc < 0 || l2.pos != 0) break;
    }

    release_lock (&cckd->filelock);
    return (l1x * 256 + l2x + dev->ckdheads) / dev->ckdheads;
}

/*-------------------------------------------------------------------*/
/* Compressed fba read block(s)                                      */
/*-------------------------------------------------------------------*/
int cfba_read_block (DEVBLK *dev, int blkgrp, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             cache;                  /* New active cache entry    */
BYTE           *cbuf;                   /* -> cache buffer           */
BYTE           *newbuf;                 /* Uncompressed buffer       */
int             len;                    /* Compressed length         */
int             maxlen;                 /* Size for cache entry      */

    cckd = dev->cckd_ext;

    if (dev->cache >= 0)
        cbuf = cache_getbuf (CACHE_DEVBUF, dev->cache, 0);
    else
        cbuf = NULL;
    maxlen = CFBA_BLOCK_SIZE + CKDDASD_TRKHDR_SIZE;

    /* Return if reading the same track image */
    if (blkgrp == dev->bufcur && dev->cache >= 0)
    {
        /* Block group image may be compressed */
        if ((cbuf[0] & CCKD_COMPRESS_MASK) != 0
         && (cbuf[0] & dev->comps) == 0)
        {
#if 0
            /* Return if synchronous i/o */
            if (dev->syncio_active)
            {
                cckd_trace (dev, "read blkgrp  %d syncio compressed\n",
                            blkgrp);
                cckdblk.stats_synciomisses++;
                dev->syncio_retry = 1;
                return -1;
            }
#endif
            len = cache_getval(CACHE_DEVBUF, dev->cache) + CKDDASD_TRKHDR_SIZE;
            newbuf = cckd_uncompress (dev, cbuf, len, maxlen, blkgrp);
            if (newbuf == NULL) {
                dev->sense[0] = SENSE_EC;
                *unitstat = CSW_CE | CSW_DE | CSW_UC;
                dev->bufcur = dev->cache = -1;
                return -1;
            }
            cache_setbuf (CACHE_DEVBUF, dev->cache, newbuf, maxlen);
            cbuf = newbuf;
            dev->buf     = newbuf + CKDDASD_TRKHDR_SIZE;
            dev->buflen  = CFBA_BLOCK_SIZE;
            cache_setval (CACHE_DEVBUF, dev->cache, dev->buflen);
            dev->bufsize = cache_getlen (CACHE_DEVBUF, dev->cache);
            dev->bufupd  = 0;
            cckd_trace (dev, "read bkgrp  %d uncompressed len %d\n",
                        blkgrp, dev->buflen);
        }

        dev->comp = cbuf[0] & CCKD_COMPRESS_MASK;

        return 0;
    }

    cckd_trace (dev, "read blkgrp  %d (%s)\n", blkgrp,
                dev->syncio_active ? "synchronous" : "asynchronous");

    /* Read the new blkgrp */
    dev->bufupd = 0;
    *unitstat = 0;
    cache = cckd_read_trk (dev, blkgrp, 0, unitstat);
    if (cache < 0)
    {
        dev->bufcur = dev->cache = -1;
        return -1;
    }
    dev->cache    = cache;
    cbuf          = cache_getbuf (CACHE_DEVBUF, dev->cache, 0);
    dev->buf      = cbuf + CKDDASD_TRKHDR_SIZE;
    dev->bufcur   = blkgrp;
    dev->bufoff   = 0;
    dev->bufoffhi = CFBA_BLOCK_SIZE;
    dev->buflen   = CFBA_BLOCK_SIZE;
    cache_setval  (CACHE_DEVBUF, dev->cache, dev->buflen);
    dev->bufsize  = cache_getlen (CACHE_DEVBUF, dev->cache);
    dev->comp     = cbuf[0] & CCKD_COMPRESS_MASK;

    /* If the image is compressed then call ourself recursively
       to cause the image to get uncompressed.  This is because
      `bufcur' will match blkgrp and `comps' won't match `comp' */
    if (dev->comp != 0 && (dev->comp & dev->comps) == 0)
        rc = cfba_read_block (dev, blkgrp, unitstat);
    else
        rc = 0;

    return rc;
} /* end function cfba_read_block */

/*-------------------------------------------------------------------*/
/* Compressed fba write block(s)                                     */
/*-------------------------------------------------------------------*/
int cfba_write_block (DEVBLK *dev, int blkgrp, int off,
                      BYTE *buf, int len, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
BYTE           *cbuf;                   /* -> cache buffer           */

    cckd = dev->cckd_ext;

    if (dev->cache >= 0)
        cbuf = cache_getbuf (CACHE_DEVBUF, dev->cache, 0);
    else
        cbuf = NULL;

    /* Read the block group if it's not current or compressed.
       `dev->comps' is set to zero forcing the read routine to
       uncompress the image.                                   */
    if (blkgrp != dev->bufcur || (cbuf[0] & CCKD_COMPRESS_MASK) != 0)
    {
        dev->comps = 0;
        rc = (dev->hnd->read) (dev, blkgrp, unitstat);
        if (rc < 0)
        {
            dev->bufcur = dev->cache = -1;
            return -1;
        }
    }

    /* Copy the data into the buffer */
    if (buf) memcpy (dev->buf + off, buf, len);

    /* Update the cache entry */
    cache_setflag (CACHE_DEVBUF, dev->cache, ~0, CCKD_CACHE_UPDATED|CCKD_CACHE_USED);
    cckd->updated = 1;

    /* Notify the shared server of the update */
    if (!dev->bufupd)
    {
        dev->bufupd = 1;
        shared_update_notify (dev, blkgrp);
    }

    return len;

} /* end function cfba_write_block */

/*-------------------------------------------------------------------*/
/* Return used blocks                                                */
/*-------------------------------------------------------------------*/
int cfba_used (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             l1x, l2x;               /* Lookup table indexes      */
int             sfx;                    /* Shadow file suffix        */
CCKD_L2ENT      l2;                     /* Copied level 2 entry      */

    cckd = dev->cckd_ext;
    obtain_lock (&cckd->filelock);

    /* Find the last used level 1 table entry */
    for (l1x = cckd->cdevhdr[0].numl1tab - 1; l1x > 0; l1x--)
    {
        sfx = cckd->sfn;
        while (cckd->l1[sfx][l1x] == 0xffffffff && sfx > 0) sfx--;
        if (cckd->l1[sfx][l1x]) break;
    }

    /* Find the last used level 2 table entry */
    for (l2x = 255; l2x >= 0; l2x--)
    {
        rc = cckd_read_l2ent (dev, &l2, l1x * 256 + l2x);
        if (rc < 0 || l2.pos != 0) break;
    }

    release_lock (&cckd->filelock);
    return (l1x * 256 + l2x + CFBA_BLOCK_NUM) / CFBA_BLOCK_NUM;
}

/*-------------------------------------------------------------------*/
/* Read a track image                                                */
/*                                                                   */
/* This routine can be called by the i/o thread (`ra' == 0) or       */
/* by readahead threads (0 < `ra' <= cckdblk.ramax).                 */
/*                                                                   */
/*-------------------------------------------------------------------*/
int cckd_read_trk(DEVBLK *dev, int trk, int ra, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             fnd;                    /* Cache index for hit       */
int             lru;                    /* Oldest unused cache index */
int             len;                    /* Length of track image     */
int             maxlen;                 /* Length for buffer         */
int             curtrk = -1;            /* Current track (at entry)  */
U16             devnum;                 /* Device number             */
U32             oldtrk;                 /* Stolen track number       */
U32             flag;                   /* Cache flag                */
BYTE           *buf;                    /* Read buffer               */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "%d rdtrk     %d\n", ra, trk);

    maxlen = cckd->ckddasd ? dev->ckdtrksz
                           : CFBA_BLOCK_SIZE + CKDDASD_TRKHDR_SIZE;

    if (!ra) obtain_lock (&cckd->iolock);

    cache_lock (CACHE_DEVBUF);

    /* Inactivate the old entry */
    if (!ra)
    {
        curtrk = dev->bufcur;
        if (dev->cache >= 0)
            cache_setflag(CACHE_DEVBUF, dev->cache, ~CCKD_CACHE_ACTIVE, 0);
        dev->bufcur = dev->cache = -1;
    }

cckd_read_trk_retry:

    /* scan the cache array for the track */
    fnd = cache_lookup (CACHE_DEVBUF, CCKD_CACHE_SETKEY(dev->devnum, trk), &lru);

    /* check for cache hit */
    if (fnd >= 0)
    {
        if (ra) /* readahead doesn't care about a cache hit */
        {   cache_unlock (CACHE_DEVBUF);
            return fnd;
        }

        /* If synchronous I/O and I/O is active then return
           with syncio_retry bit on */
        if (dev->syncio_active)
        {
            if (cache_getflag(CACHE_DEVBUF, fnd) & CCKD_CACHE_IOBUSY)
            {
                cckd_trace (dev, "%d rdtrk[%d] %d syncio %s\n", ra, fnd, trk,
                            cache_getflag(CACHE_DEVBUF, fnd) & CCKD_CACHE_READING ?
                            "reading" : "writing");
                cckdblk.stats_synciomisses++;
                dev->syncio_retry = 1;
                cache_unlock (CACHE_DEVBUF);
                release_lock (&cckd->iolock);
                return -1;
            }
            else cckdblk.stats_syncios++;
        }

        /* Mark the new entry active */
        cache_setflag(CACHE_DEVBUF, fnd, ~0, CCKD_CACHE_ACTIVE | CCKD_CACHE_USED);
        cache_setage(CACHE_DEVBUF, fnd);

        /* If the entry is pending write then change it to `updated' */
        if (cache_getflag(CACHE_DEVBUF, fnd) & CCKD_CACHE_WRITE)
        {
            cache_setflag(CACHE_DEVBUF, fnd, ~CCKD_CACHE_WRITE, CCKD_CACHE_UPDATED);
            cckd->wrpending--;
            if (cckd->iowaiters && !cckd->wrpending)
                broadcast_condition (&cckd->iocond);
        }
        buf = cache_getbuf(CACHE_DEVBUF, fnd, 0);

        cache_unlock (CACHE_DEVBUF);

        cckd_trace (dev, "%d rdtrk[%d] %d cache hit buf %p:%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                    ra, fnd, trk, buf, buf[0], buf[1], buf[2], buf[3], buf[4]);

        cckdblk.stats_switches++;  cckd->switches++;
        cckdblk.stats_cachehits++; cckd->cachehits++;

        /* if read/write is in progress then wait for it to finish */
        while (cache_getflag(CACHE_DEVBUF, fnd) & CCKD_CACHE_IOBUSY)
        {
            cckdblk.stats_iowaits++;
            cckd_trace (dev, "%d rdtrk[%d] %d waiting for %s\n", ra, fnd, trk,
                        cache_getflag(CACHE_DEVBUF, fnd) & CCKD_CACHE_READING ?
                        "read" : "write");
            cache_setflag (CACHE_DEVBUF, fnd, ~0, CCKD_CACHE_IOWAIT);
            cckd->iowaiters++;
            wait_condition (&cckd->iocond, &cckd->iolock);
            cckd->iowaiters--;
            cache_setflag (CACHE_DEVBUF, fnd, ~CCKD_CACHE_IOWAIT, 0);
            cckd_trace (dev, "%d rdtrk[%d] %d io wait complete\n",
                        ra, fnd, trk);
        }

        release_lock (&cckd->iolock);

        /* Asynchrously schedule readaheads */
        if (curtrk > 0 && trk > curtrk && trk <= curtrk + 2)
            cckd_readahead (dev, trk);

        return fnd;

    } /* cache hit */

    /* If not readahead and synchronous I/O then retry */
    if (!ra && dev->syncio_active)
    {
        cache_unlock(CACHE_DEVBUF);
        release_lock (&cckd->iolock);
        cckd_trace (dev, "%d rdtrk[%d] %d syncio cache miss\n", ra, lru, trk);
        cckdblk.stats_synciomisses++;
        dev->syncio_retry = 1;
        return -1;
    }

    cckd_trace (dev, "%d rdtrk[%d] %d cache miss\n", ra, lru, trk);

    /* If no cache entry was stolen, then flush all outstanding writes.
       This requires us to release our locks.  cache_wait should be
       called with only the cache_lock held.  Fortunately, cache waits
       occur very rarely. */
    if (lru < 0) /* No available entry to be stolen */
    {
        cckd_trace (dev, "%d rdtrk[%d] %d no available cache entry\n",
                    ra, lru, trk);
        cache_unlock (CACHE_DEVBUF);
        if (!ra) release_lock (&cckd->iolock);
        cckd_flush_cache_all();
        cache_lock (CACHE_DEVBUF);
        cckdblk.stats_cachewaits++;
        cache_wait (CACHE_DEVBUF);
        if (!ra)
        {
            cache_unlock (CACHE_DEVBUF);
            obtain_lock (&cckd->iolock);
            cache_lock (CACHE_DEVBUF);
        }
        goto cckd_read_trk_retry;
    }

    CCKD_CACHE_GETKEY(lru, devnum, oldtrk);
    if (devnum != 0)
    {
        cckd_trace (dev, "%d rdtrk[%d] %d dropping %4.4X:%d from cache\n",
                    ra, lru, trk, devnum, oldtrk);
        if (!(cache_getflag(CACHE_DEVBUF, lru) & CCKD_CACHE_USED))
        {
            cckdblk.stats_readaheadmisses++;  cckd->misses++;
        }
    }

    /* Initialize the entry */
    cache_setkey(CACHE_DEVBUF, lru, CCKD_CACHE_SETKEY(dev->devnum, trk));
    cache_setflag(CACHE_DEVBUF, lru, 0, CCKD_CACHE_READING);
    cache_setage(CACHE_DEVBUF, lru);
    cache_setval(CACHE_DEVBUF, lru, 0);
    if (!ra)
    {
        cckdblk.stats_switches++; cckd->switches++;
        cckdblk.stats_cachemisses++;
        cache_setflag(CACHE_DEVBUF, lru, ~0, CCKD_CACHE_ACTIVE|CCKD_CACHE_USED);
    }
    cache_setflag(CACHE_DEVBUF, lru, ~CACHE_TYPE,
                  cckd->ckddasd ? DEVBUF_TYPE_CCKD : DEVBUF_TYPE_CFBA);
    buf = cache_getbuf(CACHE_DEVBUF, lru, maxlen);

    cckd_trace (dev, "%d rdtrk[%d] %d buf %p len %d\n",
                ra, lru, trk, buf, cache_getlen(CACHE_DEVBUF, lru));

    cache_unlock (CACHE_DEVBUF);

    if (!ra) release_lock (&cckd->iolock);

    /* Asynchronously schedule readaheads */
    if (!ra && curtrk > 0 && trk > curtrk && trk <= curtrk + 2)
        cckd_readahead (dev, trk);

    /* Clear the buffer if batch mode */
    if (dev->batch) memset(buf, 0, maxlen);

    /* Read the track image */
    obtain_lock (&cckd->filelock);
    len = cckd_read_trkimg (dev, buf, trk, unitstat);
    release_lock (&cckd->filelock);
    cache_setval (CACHE_DEVBUF, lru, len);

    obtain_lock (&cckd->iolock);

    /* Turn off the READING bit */
    cache_lock (CACHE_DEVBUF);
    flag = cache_setflag(CACHE_DEVBUF, lru, ~CCKD_CACHE_READING, 0);
    cache_unlock (CACHE_DEVBUF);

    /* Wakeup other thread waiting for this read */
    if (cckd->iowaiters && (flag & CCKD_CACHE_IOWAIT))
    {   cckd_trace (dev, "%d rdtrk[%d] %d signalling read complete\n",
                    ra, lru, trk);
        broadcast_condition (&cckd->iocond);
    }

    release_lock (&cckd->iolock);

    if (ra)
    {
        cckdblk.stats_readaheads++; cckd->readaheads++;
    }

    cckd_trace (dev, "%d rdtrk[%d] %d complete buf %p:%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                ra, lru, trk, buf, buf[0], buf[1], buf[2], buf[3], buf[4]);

    if (cache_busy_percent(CACHE_DEVBUF) > 80) cckd_flush_cache_all();

    return lru;

} /* end function cckd_read_trk */

/*-------------------------------------------------------------------*/
/* Schedule asynchronous readaheads                                  */
/*-------------------------------------------------------------------*/
void cckd_readahead (DEVBLK *dev, int trk)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             i, r;                   /* Indexes                   */
TID             tid;                    /* Readahead thread id       */

    cckd = dev->cckd_ext;

    if (cckdblk.ramax < 1 || cckdblk.readaheads < 1)
        return;

    obtain_lock (&cckdblk.ralock);

    /* Scan the cache to see if the tracks are already there */
    memset(cckd->ralkup, 0, sizeof(cckd->ralkup));
    cckd->ratrk = trk;
    cache_lock(CACHE_DEVBUF);
    cache_scan(CACHE_DEVBUF, cckd_readahead_scan, dev);
    cache_unlock(CACHE_DEVBUF);

    /* Scan the queue to see if the tracks are already there */
    for (r = cckdblk.ra1st; r >= 0; r = cckdblk.ra[r].next)
        if (cckdblk.ra[r].dev == dev)
        {
            i = cckdblk.ra[r].trk - trk;
            if (i > 0 && i <= cckdblk.readaheads)
                cckd->ralkup[i-1] = 1;
        }

    /* Queue the tracks to the readahead queue */
    for (i = 1; i <= cckdblk.readaheads && cckdblk.rafree >= 0; i++)
    {
        if (cckd->ralkup[i-1]) continue;
        if (trk + i >= dev->ckdtrks) break;
        r = cckdblk.rafree;
        cckdblk.rafree = cckdblk.ra[r].next;
        if (cckdblk.ralast < 0)
        {
            cckdblk.ra1st = cckdblk.ralast = r;
            cckdblk.ra[r].prev = cckdblk.ra[r].next = -1;
        }
        else
        {
            cckdblk.ra[cckdblk.ralast].next = r;
            cckdblk.ra[r].prev = cckdblk.ralast;
            cckdblk.ra[r].next = -1;
            cckdblk.ralast = r;
        }
        cckdblk.ra[r].trk = trk + i;
        cckdblk.ra[r].dev = dev;
    }

    /* Schedule the readahead if any are pending */
    if (cckdblk.ra1st >= 0)
    {
        if (cckdblk.rawaiting)
            signal_condition (&cckdblk.racond);
        else if (cckdblk.ras < cckdblk.ramax)
            create_thread (&tid, JOINABLE, cckd_ra, NULL, "cckd_ra");
    }

    release_lock (&cckdblk.ralock);

} /* end function cckd_readahead */

int cckd_readahead_scan (int *answer, int ix, int i, void *data)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
U16             devnum;                 /* Cached device number      */
U32             trk;                    /* Cached track              */
DEVBLK         *dev = data;             /* -> device block           */
int             k;                      /* Index                     */

    UNREFERENCED(answer);
    UNREFERENCED(ix);
    cckd = dev->cckd_ext;
    CCKD_CACHE_GETKEY(i, devnum, trk);
    if (devnum == dev->devnum)
    {
        k = (int)trk - cckd->ratrk;
        if (k > 0 && k <= cckdblk.readaheads)
            cckd->ralkup[k-1] = 1;
    }
    return 0;
}

/*-------------------------------------------------------------------*/
/* Asynchronous readahead thread                                     */
/*-------------------------------------------------------------------*/
void *cckd_ra ()
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
DEVBLK         *dev;                    /* Readahead devblk          */
int             trk;                    /* Readahead track           */
int             ra;                     /* Readahead index           */
int             r;                      /* Readahead queue index     */
TID             tid;                    /* Readahead thread id       */

    obtain_lock (&cckdblk.ralock);
    ra = ++cckdblk.ras;

    /* Return without messages if too many already started */
    if (ra > cckdblk.ramax)
    {
        --cckdblk.ras;
        release_lock (&cckdblk.ralock);
        return NULL;
    }

    if (!cckdblk.batch)
    {
        logmsg (_("HHCCD001I Readahead thread %d started: tid="TIDPAT", pid=%d\n"),
            ra, thread_id(), getpid());
    }

    while (ra <= cckdblk.ramax)
    {
        if (cckdblk.ra1st < 0)
        {
            cckdblk.rawaiting++;
            wait_condition (&cckdblk.racond, &cckdblk.ralock);
            cckdblk.rawaiting--;
        }

        /* Possibly shutting down if no writes pending */
        if (cckdblk.ra1st < 0) continue;

        r = cckdblk.ra1st;
        trk = cckdblk.ra[r].trk;
        dev = cckdblk.ra[r].dev;
        cckd = dev->cckd_ext;

        /* Requeue the 1st entry to the readahead free queue */
        cckdblk.ra1st = cckdblk.ra[r].next;
        if (cckdblk.ra[r].next > -1)
            cckdblk.ra[cckdblk.ra[r].next].prev = -1;
        else cckdblk.ralast = -1;
        cckdblk.ra[r].next = cckdblk.rafree;
        cckdblk.rafree = r;

        /* Schedule the other readaheads if any are still pending */
        if (cckdblk.ra1st)
        {
            if (cckdblk.rawaiting)
                signal_condition (&cckdblk.racond);
            else if (cckdblk.ras < cckdblk.ramax)
                create_thread (&tid, JOINABLE, cckd_ra, dev, "cckd_ra");
        }

        if (!cckd || cckd->stopping || cckd->merging) continue;

        cckd->ras++;
        release_lock (&cckdblk.ralock);

        /* Read the readahead track */
        cckd_read_trk (dev, trk, ra, NULL);

        obtain_lock (&cckdblk.ralock);
        cckd->ras--;
    }

    if (!cckdblk.batch)
    logmsg (_("HHCCD011I Readahead thread %d stopping: tid="TIDPAT", pid=%d\n"),
            ra, thread_id(), getpid());
    --cckdblk.ras;
    if (!cckdblk.ras) signal_condition(&cckdblk.termcond);
    release_lock(&cckdblk.ralock);
    return NULL;

} /* end thread cckd_ra_thread */

/*-------------------------------------------------------------------*/
/* Flush updated cache entries for a device                          */
/*                                                                   */
/* Caller holds the cckd->iolock                                     */
/* cckdblk.wrlock then cache_lock is obtained and released           */
/*-------------------------------------------------------------------*/
void cckd_flush_cache(DEVBLK *dev)
{
int             rc;                     /* Return code               */
TID             tid;                    /* Writer thread id          */

    /* Scan cache for updated cache entries */
    obtain_lock (&cckdblk.wrlock);
    cache_lock (CACHE_DEVBUF);
    rc = cache_scan (CACHE_DEVBUF, cckd_flush_cache_scan, dev);
    cache_unlock (CACHE_DEVBUF);

    /* Schedule the writer if any writes are pending */
    if (cckdblk.wrpending)
    {
        if (cckdblk.wrwaiting)
            signal_condition (&cckdblk.wrcond);
        else if (cckdblk.wrs < cckdblk.wrmax)
        {
            create_thread (&tid, JOINABLE, cckd_writer, NULL, "cckd_writer");
        }
    }
    release_lock (&cckdblk.wrlock);
}
int cckd_flush_cache_scan (int *answer, int ix, int i, void *data)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
U16             devnum;                 /* Cached device number      */
U32             trk;                    /* Cached track              */
DEVBLK         *dev = data;             /* -> device block           */

    UNREFERENCED(answer);
    cckd = dev->cckd_ext;
    CCKD_CACHE_GETKEY(i, devnum, trk);
    if ((cache_getflag(ix,i) & CACHE_BUSY) == CCKD_CACHE_UPDATED
     && dev->devnum == devnum)
    {
        cache_setflag (ix, i, ~CCKD_CACHE_UPDATED, CCKD_CACHE_WRITE);
        ++cckd->wrpending;
        ++cckdblk.wrpending;
        cckd_trace (dev, "flush file[%d] cache[%d] %4.4X trk %d\n",
                    cckd->sfn, i, devnum, trk);
    }
    return 0;
}
void cckd_flush_cache_all()
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
DEVBLK         *dev = NULL;             /* -> device block           */

    cckd_lock_devchain(0);
    for (dev = cckdblk.dev1st; dev; dev = cckd->devnext)
    {
        cckd = dev->cckd_ext;
        obtain_lock (&cckd->iolock);
        if (!cckd->merging && !cckd->stopping)
            cckd_flush_cache(dev);
        release_lock (&cckd->iolock);
    }
    cckd_unlock_devchain();
}

/*-------------------------------------------------------------------*/
/* Purge cache entries for a device                                  */
/*                                                                   */
/* Caller holds the iolock                                           */
/* cache_lock is obtained and released                               */
/*-------------------------------------------------------------------*/
void cckd_purge_cache(DEVBLK *dev)
{

    /* Scan cache and purge entries */
    cache_lock (CACHE_DEVBUF);
    cache_scan (CACHE_DEVBUF, cckd_purge_cache_scan, dev);
    cache_unlock (CACHE_DEVBUF);
}
int cckd_purge_cache_scan (int *answer, int ix, int i, void *data)
{
U16             devnum;                 /* Cached device number      */
U32             trk;                    /* Cached track              */
DEVBLK         *dev = data;             /* -> device block           */

    UNREFERENCED(answer);
    CCKD_CACHE_GETKEY(i, devnum, trk);
    if (dev->devnum == devnum)
    {
        cache_release (ix, i, 0);
        cckd_trace (dev, "purge cache[%d] %4.4X trk %d purged\n",
                    i, devnum, trk);
    }
    return 0;
}

/*-------------------------------------------------------------------*/
/* Writer thread                                                     */
/*-------------------------------------------------------------------*/
void *cckd_writer(void *arg)
{
DEVBLK         *dev;                    /* Device block              */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             writer;                 /* Writer identifier         */
int             o;                      /* Cache entry found         */
U16             devnum;                 /* Device number             */
BYTE           *buf;                    /* Buffer                    */
BYTE           *bufp;                   /* Buffer to be written      */
int             len, bufl;              /* Buffer lengths            */
int             trk;                    /* Track number              */
int             comp;                   /* Compression algorithm     */
int             parm;                   /* Compression parameter     */
TID             tid;                    /* Writer thead id           */
U32             flag;                   /* Cache flag                */
static char    *compress[] = {"none", "zlib", "bzip2"};
BYTE            buf2[65536];            /* Compress buffer           */

    UNREFERENCED(arg);

#ifndef WIN32
    /* Set writer priority just below cpu priority to mimimize the
       compression effect */
    if(cckdblk.wrprio >= 0)
        setpriority (PRIO_PROCESS, 0, cckdblk.wrprio);
#endif

    obtain_lock (&cckdblk.wrlock);

    writer = ++cckdblk.wrs;

    /* Return without messages if too many already started */
    if (writer > cckdblk.wrmax)
    {
        --cckdblk.wrs;
        release_lock (&cckdblk.wrlock);
        return NULL;
    }

    if (!cckdblk.batch)
    {
        logmsg (_("HHCCD002I Writer thread %d started: tid="TIDPAT", pid=%d\n"),
            writer, thread_id(), getpid());
    }

    while (writer <= cckdblk.wrmax || cckdblk.wrpending)
    {
        /* Wait for work */
        if (cckdblk.wrpending == 0)
        {
            cckdblk.wrwaiting++;
            wait_condition (&cckdblk.wrcond, &cckdblk.wrlock);
            cckdblk.wrwaiting--;
        }

        /* Scan the cache for the oldest pending write */
        cache_lock (CACHE_DEVBUF);
        o = cache_scan (CACHE_DEVBUF, cckd_writer_scan, NULL);

        /* Possibly shutting down if no writes pending */
        if (o < 0)
        {
            cache_unlock (CACHE_DEVBUF);
            cckdblk.wrpending = 0;
            continue;
        }
        cache_setflag (CACHE_DEVBUF, o, ~CCKD_CACHE_WRITE, CCKD_CACHE_WRITING);
        cache_unlock (CACHE_DEVBUF);

        /* Schedule the other writers if any writes are still pending */
        cckdblk.wrpending--;
        if (cckdblk.wrpending)
        {
            if (cckdblk.wrwaiting)
                signal_condition (&cckdblk.wrcond);
            else if (cckdblk.wrs < cckdblk.wrmax)
            {
                create_thread (&tid, JOINABLE, cckd_writer, NULL, "cckd_writer");
            }
        }
        release_lock (&cckdblk.wrlock);

        /* Prepare to compress */
        CCKD_CACHE_GETKEY(o, devnum, trk);
        dev = cckd_find_device_by_devnum (devnum);
        cckd = dev->cckd_ext;
        buf = cache_getbuf(CACHE_DEVBUF, o, 0);
        len = cckd_trklen (dev, buf);
        comp = len < CCKD_COMPRESS_MIN ? CCKD_COMPRESS_NONE
             : cckdblk.comp == 0xff ? cckd->cdevhdr[cckd->sfn].compress
             : cckdblk.comp;
        parm = cckdblk.compparm < 0
             ? cckd->cdevhdr[cckd->sfn].compress_parm
             : cckdblk.compparm;

        cckd_trace (dev, "%d wrtrk[%d] %d len %d buf %p:%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                    writer, o, trk, len, buf, buf[0], buf[1],buf[2],buf[3],buf[4]);

        /* Compress the image if not null */
        if ((len = cckd_check_null_trk (dev, buf, trk, len)) > CKDDASD_NULLTRK_FMTMAX)
        {
            /* Stress adjustments */
            if ((cache_waiters(CACHE_DEVBUF) || cache_busy(CACHE_DEVBUF) > 90)
             && !cckdblk.nostress)
            {
                cckdblk.stats_stresswrites++;
                comp = len < CCKD_STRESS_MINLEN ?
                       CCKD_COMPRESS_NONE : CCKD_STRESS_COMP;
                parm = cache_busy(CACHE_DEVBUF) <= 95 ?
                       CCKD_STRESS_PARM1 : CCKD_STRESS_PARM2;
            }

            /* Compress the track image */
            cckd_trace (dev, "%d wrtrk[%d] %d comp %s parm %d\n",
                        writer, o, trk, compress[comp], parm);
            bufp = (BYTE *)&buf2;
            bufl = cckd_compress(dev, &bufp, buf, len, comp, parm);
            cckd_trace (dev, "%d wrtrk[%d] %d compressed length %d\n",
                        writer, o, trk, bufl);
        }
        else
        {
            bufp = buf;
            bufl = len;
        }

        obtain_lock (&cckd->filelock);

        /* Turn on read-write header bits if not already on */
        if (!(cckd->cdevhdr[cckd->sfn].options & CCKD_OPENED))
        {
            cckd->cdevhdr[cckd->sfn].options |= (CCKD_OPENED | CCKD_ORDWR);
            cckd_write_chdr (dev);
        }

        /* Write the track image */
        cckd_write_trkimg (dev, bufp, bufl, trk, CCKD_SIZE_ANY);

        release_lock (&cckd->filelock);

        /* Schedule the garbage collector */
        if (cckdblk.gcs < cckdblk.gcmax)
            create_thread (&tid, JOINABLE, cckd_gcol, NULL, "cckd_gcol");

        obtain_lock (&cckd->iolock);
        cache_lock (CACHE_DEVBUF);
        flag = cache_setflag (CACHE_DEVBUF, o, ~CCKD_CACHE_WRITING, 0);
        cache_unlock (CACHE_DEVBUF);
        cckd->wrpending--;
        if (cckd->iowaiters && ((flag & CCKD_CACHE_IOWAIT) || !cckd->wrpending))
        {   cckd_trace (dev, "writer[%d] cache[%2.2d] %d signalling write complete\n",
                        writer, o, trk);
            broadcast_condition (&cckd->iocond);
        }
        release_lock(&cckd->iolock);

        cckd_trace (dev, "%d wrtrk[%2.2d] %d complete flags:%8.8x\n",
                    writer, o, trk, cache_getflag(CACHE_DEVBUF,o));

        obtain_lock(&cckdblk.wrlock);
    }

    if (!cckdblk.batch)
    logmsg (_("HHCCD012I Writer thread %d stopping: tid="TIDPAT", pid=%d\n"),
            writer, thread_id(), getpid());
    cckdblk.wrs--;
    if (cckdblk.wrs == 0) signal_condition(&cckdblk.termcond);
    release_lock(&cckdblk.wrlock);
    return NULL;
} /* end thread cckd_writer */

int cckd_writer_scan (int *o, int ix, int i, void *data)
{
    UNREFERENCED(data);
    if ((cache_getflag(ix,i) & DEVBUF_TYPE_COMP)
     && (cache_getflag(ix,i) & CCKD_CACHE_WRITE)
     && (*o == -1 || cache_getage(ix, i) < cache_getage(ix, *o)))
        *o = i;
    return 0;
}

/*-------------------------------------------------------------------*/
/* Debug routine for checking the free space array                   */
/*-------------------------------------------------------------------*/

void cckd_chk_space(DEVBLK *dev)
{
#if 1
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* Shadow file index         */
int             err = 0, n = 0, i, p;
size_t  largest=0;
size_t  total=0;
off_t           fpos;

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    p = -1;
    fpos = cckd->cdevhdr[sfx].free;
    for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
    {
        n++; total += cckd->free[i].len;
        if (n > cckd->freenbr) break;
        if (cckd->free[i].prev != p)
            err = 1;
        if (cckd->free[i].next >= 0)
        {
            if (fpos + cckd->free[i].len > cckd->free[i].pos)
                err = 1;
        }
        else
        {
            if (fpos + cckd->free[i].len > cckd->cdevhdr[sfx].size)
                err = 1;
        }
        if (cckd->free[i].pending == 0 && cckd->free[i].len > largest)
            largest = cckd->free[i].len;
        fpos = cckd->free[i].pos;
        p = i;
    }

    if (err
     || (cckd->cdevhdr[sfx].free != 0 && cckd->cdevhdr[sfx].free_number == 0)
     || (cckd->cdevhdr[sfx].free == 0 && cckd->cdevhdr[sfx].free_number != 0)
     || (n != cckd->cdevhdr[sfx].free_number)
     || (total != cckd->cdevhdr[sfx].free_total - cckd->cdevhdr[sfx].free_imbed)
     || (cckd->freelast != p)
     || (largest != cckd->cdevhdr[sfx].free_largest)
       )
    {
        cckd_trace (dev, "cdevhdr[%d] size   %10d used   %10d free   0x%8.8x\n",
                    sfx,cckd->cdevhdr[sfx].size,cckd->cdevhdr[sfx].used,
                    cckd->cdevhdr[sfx].free);
        cckd_trace (dev, "           nbr   %10d total  %10d imbed  %10d largest %10d\n",
                    cckd->cdevhdr[sfx].free_number,
                    cckd->cdevhdr[sfx].free_total,cckd->cdevhdr[sfx].free_imbed,
                    cckd->cdevhdr[sfx].free_largest);
        cckd_trace (dev, "free %p nbr %d 1st %d last %d avail %d\n",
                    cckd->free,cckd->freenbr,cckd->free1st,
                    cckd->freelast,cckd->freeavail);
        cckd_trace (dev, "found nbr %d total %ld largest %ld\n",n,(long)total,(long)largest);
        fpos = cckd->cdevhdr[sfx].free;
        for (n = 0, i = cckd->free1st; i >= 0; i = cckd->free[i].next)
        {
            if (++n > cckd->freenbr) break;
            cckd_trace (dev, "%4d: [%4d] prev[%4d] next[%4d] pos %8.8" I64_FMT "x len %8d %8.8" I64_FMT "x pend %d\n",
                        n, i, cckd->free[i].prev, cckd->free[i].next,
                        (long long)fpos, cckd->free[i].len,
                        (long long)fpos + cckd->free[i].len, cckd->free[i].pending);
            fpos = cckd->free[i].pos;
        }
        cckd_print_itrace();
    }
#endif
} /* end function cckd_chk_space */

/*-------------------------------------------------------------------*/
/* Get file space                                                    */
/*-------------------------------------------------------------------*/
off_t cckd_get_space(DEVBLK *dev, int *size, int flags)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             i,p,n;                  /* Free space indexes        */
int             len2;                   /* Other lengths             */
off_t           fpos;                   /* Free space offset         */
unsigned int    flen;                   /* Free space size           */
int             sfx;                    /* Shadow file index         */
int             len;                    /* Requested length          */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    len = *size;

    if (flags & CCKD_L2SPACE)
    {
        flags |= CCKD_SIZE_EXACT;
        len = *size = CCKD_L2TAB_SIZE;
    }

    cckd_trace (dev, "get_space len %d largest %d flags 0x%2.2x\n",
                len, cckd->cdevhdr[sfx].free_largest, flags);

    if (len <= CKDDASD_NULLTRK_FMTMAX)
        return 0;

    if (!cckd->free)
        cckd_read_fsp (dev);

    len2 = len + CCKD_FREEBLK_SIZE;

    /* Get space at the end if no space is large enough */
    if (len2 > (int)cckd->cdevhdr[sfx].free_largest
     && len != (int)cckd->cdevhdr[sfx].free_largest)
    {

cckd_get_space_atend:

        fpos = (off_t)cckd->cdevhdr[sfx].size;
        if ((long long)(fpos + len) > cckd->maxsize)
        {
            logmsg (_("HHCCD102E %4.4X file[%d] get space error, size exceeds %lldM\n"),
                    dev->devnum, sfx, (cckd->maxsize >> 20) + 1);
            return -1;
        }
        cckd->cdevhdr[sfx].size += len;
        cckd->cdevhdr[sfx].used += len;

        cckd_trace (dev, "get_space atend 0x%" I64_FMT "x len %d\n",(long long)fpos, len);

        return fpos;
    }

    /* Scan free space chain */
    fpos = (off_t)cckd->cdevhdr[sfx].free;
    for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
    {
        if (cckd->free[i].pending == 0
         && (len2 <= (int)cckd->free[i].len || len == (int)cckd->free[i].len)
         && ((flags & CCKD_L2SPACE) || fpos >= cckd->l2bounds))
            break;
        fpos = (off_t)cckd->free[i].pos;
    }

    /* This can happen if largest comes before l2bounds */
    if (i < 0) goto cckd_get_space_atend; 

    flen = cckd->free[i].len;
    p = cckd->free[i].prev;
    n = cckd->free[i].next;

    /*
     * If `CCKD_SIZE_ANY' bit is set and the left over space is small
     * enough, then use the entire free space
     */
    if ((flags & CCKD_SIZE_ANY) && flen <= cckd->freemin)
        *size = (int)flen;

    /* Remove the new space from free space */
    if (*size < (int)flen)
    {
        cckd->free[i].len -= *size;
        if (p >= 0)
            cckd->free[p].pos += *size;
        else
            cckd->cdevhdr[sfx].free += *size;
    }
    else
    {
        cckd->cdevhdr[sfx].free_number--;

        /* Remove the free space entry from the chain */
        if (p >= 0)
        {
            cckd->free[p].pos = cckd->free[i].pos;
            cckd->free[p].next = n;
        }
        else
        {
            cckd->cdevhdr[sfx].free = cckd->free[i].pos;
            cckd->free1st = n;
        }

        if (n >= 0)
            cckd->free[n].prev = p;
        else
            cckd->freelast = p;

        /* Add entry to the available queue */
        cckd->free[i].next = cckd->freeavail;
        cckd->freeavail = i;
    }

    /* Find the largest free space if we got the largest */
    if (flen >= cckd->cdevhdr[sfx].free_largest)
    {
        int i;
        cckd->cdevhdr[sfx].free_largest = 0;
        for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
            if (cckd->free[i].len > cckd->cdevhdr[sfx].free_largest
             && cckd->free[i].pending == 0)
                cckd->cdevhdr[sfx].free_largest = cckd->free[i].len;
    }

    /* Update free space stats */
    cckd->cdevhdr[sfx].used += len;
    cckd->cdevhdr[sfx].free_total -= len;

    cckd->cdevhdr[sfx].free_imbed += *size - len;

    cckd_trace (dev, "get_space found 0x%" I64_FMT "x len %d size %d\n",
                (long long)fpos, len, *size);

    return fpos;

} /* end function cckd_get_space */

/*-------------------------------------------------------------------*/
/* Release file space                                                */
/*-------------------------------------------------------------------*/
void cckd_rel_space(DEVBLK *dev, off_t pos, int len, int size)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* Shadow file index         */
off_t           ppos, npos;             /* Prev/next free offsets    */
int             i, p, n;                /* Free space indexes        */
int             pending;                /* Calculated pending value  */
int             fsize = size;           /* Free space size           */

    if (len <= CKDDASD_NULLTRK_FMTMAX || pos == 0 || pos == 0xffffffff)
        return;

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "rel_space offset %" I64_FMT "x len %d size %d\n",
                (long long)pos, len, size);

    if (!cckd->free) cckd_read_fsp (dev);

//  cckd_chk_space(dev);

    /* Scan free space chain */
    ppos = -1;
    npos = cckd->cdevhdr[sfx].free;
    for (p = -1, n = cckd->free1st; n >= 0; n = cckd->free[n].next)
    {
        if (pos < npos) break;
        ppos = npos;
        npos = cckd->free[n].pos;
        p = n;
    }

    /* Calculate the `pending' value */
    pending = cckdblk.freepend >= 0 ? cckdblk.freepend : 1 + (1 - cckdblk.fsync);

    /* If possible use previous adjacent free space otherwise get an available one */
    if (p >= 0 && ppos + cckd->free[p].len == pos && cckd->free[p].pending == pending)
    {
        cckd->free[p].len += size;
        fsize = cckd->free[p].len;
    }
    else
    {
        /* Increase the size of the free space array if necessary */
        if (cckd->freeavail < 0)
        {
            cckd->freeavail = cckd->freenbr;
            cckd->freenbr += 1024;
            cckd->free = realloc ( cckd->free, cckd->freenbr * CCKD_FREEBLK_ISIZE);
            for (i = cckd->freeavail; i < cckd->freenbr; i++)
                cckd->free[i].next = i + 1;
            cckd->free[i-1].next = -1;
            cckd->freemin = CCKD_FREE_MIN_SIZE + ((cckd->freenbr >> 10) * CCKD_FREE_MIN_INCR);
        }

        /* Get an available free space entry */
        i = cckd->freeavail;
        cckd->freeavail = cckd->free[i].next;
        cckd->cdevhdr[sfx].free_number++;

        /* Update the new entry */
        cckd->free[i].prev = p;
        cckd->free[i].next = n;
        cckd->free[i].len = size;
        cckd->free[i].pending = pending;

        /* Update the previous entry */
        if (p >= 0)
        {
            cckd->free[i].pos = cckd->free[p].pos;
            cckd->free[p].pos = pos;
            cckd->free[p].next = i;
        }
        else
        {
            cckd->free[i].pos = cckd->cdevhdr[sfx].free;
            cckd->cdevhdr[sfx].free = pos;
            cckd->free1st = i;
        }

        /* Update the next entry */
        if (n >= 0)
            cckd->free[n].prev = i;
        else
            cckd->freelast = i;
    }

    /* Update the free space statistics */
    cckd->cdevhdr[sfx].used -= len;
    cckd->cdevhdr[sfx].free_total += len;
    cckd->cdevhdr[sfx].free_imbed -= size - len;
    if (!pending && (U32)fsize > cckd->cdevhdr[sfx].free_largest)
        cckd->cdevhdr[sfx].free_largest = (U32)fsize;

//  cckd_chk_space(dev);

} /* end function cckd_rel_space */

/*-------------------------------------------------------------------*/
/* Flush pending free space                                          */
/*-------------------------------------------------------------------*/
void cckd_flush_space(DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             p,i,n;                  /* Free free space indexes   */
int             sfx;                    /* Shadow file index         */
U32             ppos, pos;              /* Free space offsets        */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "flush_space nbr %d\n",cckd->cdevhdr[sfx].free_number);

    /* Make sure the free space chain is built */
    if (!cckd->free) cckd_read_fsp (dev);

//  cckd_chk_space(dev);

    if (cckd->cdevhdr[sfx].free_number == 0 || cckd->cdevhdr[sfx].free == 0)
    {
        cckd->cdevhdr[sfx].free_number = cckd->cdevhdr[sfx].free = 0;
        cckd->free1st = cckd->freelast = cckd->freeavail = -1;
    }

    pos = cckd->cdevhdr[sfx].free;
    ppos = p = -1;
    cckd->cdevhdr[sfx].free_number = cckd->cdevhdr[sfx].free_largest = 0;
    for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
    {
        /* Decrement the pending count */
        if (cckd->free[i].pending)
            --cckd->free[i].pending;

        /* Combine adjacent free spaces */
        while (pos + cckd->free[i].len == cckd->free[i].pos)
        {
            n = cckd->free[i].next;
            if (cckd->free[n].pending > cckd->free[i].pending + 1
             || cckd->free[n].pending < cckd->free[i].pending)
                break;
            cckd->free[i].pos = cckd->free[n].pos;
            cckd->free[i].len += cckd->free[n].len;
            cckd->free[i].next = cckd->free[n].next;
            cckd->free[n].next = cckd->freeavail;
            cckd->freeavail = n;
            n = cckd->free[i].next;
            if (n >= 0)
                cckd->free[n].prev = i;

        }
        ppos = pos;
        pos = cckd->free[i].pos;
        cckd->cdevhdr[sfx].free_number++;
        if (cckd->free[i].len > cckd->cdevhdr[sfx].free_largest
         && !cckd->free[i].pending)
            cckd->cdevhdr[sfx].free_largest = cckd->free[i].len;
        p = i;
    }
    cckd->freelast = p;

    cckd_trace (dev, "rel_flush_space nbr %d (after merge)\n",
                cckd->cdevhdr[sfx].free_number);

    /* If the last free space is at the end of the file then release it */
    if (p >= 0 && ppos + cckd->free[p].len == cckd->cdevhdr[sfx].size
     && !cckd->free[p].pending)
    {
        i = p;
        p = cckd->free[i].prev;

        cckd_trace (dev, "file[%d] rel_flush_space atend 0x%" I64_FMT "x len %d\n",
                    sfx, (long long)ppos, cckd->free[i].len);

        /* Remove the entry from the chain */
        if (p >= 0)
        {
            cckd->free[p].pos = 0;
            cckd->free[p].next = -1;
        }
        else
        {
            cckd->cdevhdr[sfx].free = 0;
            cckd->free1st = -1;
        }
        cckd->freelast = p;

        /* Add the entry to the available chain */
        cckd->free[i].next = cckd->freeavail;
        cckd->freeavail = i;

        /* Update the device header */
        cckd->cdevhdr[sfx].size -= cckd->free[i].len;
        cckd->cdevhdr[sfx].free_total -= cckd->free[i].len;
        cckd->cdevhdr[sfx].free_number--;
        if (cckd->free[i].len >= cckd->cdevhdr[sfx].free_largest)
        {
            cckd->cdevhdr[sfx].free_largest = 0;
            for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
                if (cckd->free[i].len > cckd->cdevhdr[sfx].free_largest
                 && cckd->free[i].pending == 0)
                    cckd->cdevhdr[sfx].free_largest = cckd->free[i].len;
        }

        /* Truncate the file */
        cckd_ftruncate (dev, sfx, (off_t)cckd->cdevhdr[sfx].size);

    } /* Release space at end of the file */

//  cckd_chk_space(dev);

} /* end function cckd_flush_space */

/*-------------------------------------------------------------------*/
/* Read compressed dasd header                                       */
/*-------------------------------------------------------------------*/
int cckd_read_chdr (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* File index                */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "file[%d] read_chdr\n", sfx);

    memset (&cckd->cdevhdr[sfx], 0, CCKDDASD_DEVHDR_SIZE);

    /* Read the device header */
    if (cckd_read (dev, sfx, CKDDASD_DEVHDR_SIZE, &cckd->cdevhdr[sfx], CCKDDASD_DEVHDR_SIZE) < 0)
        return -1;

    /* Check endian format */
    cckd->swapend[sfx] = 0;
    if ((cckd->cdevhdr[sfx].options & CCKD_BIGENDIAN) != cckd_endian())
    {
        if (cckd->open[sfx] == CCKD_OPEN_RW)
        {
            if (cckd_swapend (dev) < 0)
                return -1;
            cckd_swapend_chdr (&cckd->cdevhdr[sfx]);
        }
        else
        {
            cckd->swapend[sfx] = 1;
            cckd_swapend_chdr (&cckd->cdevhdr[sfx]);
        }
    }

    /* Set default null format */
    if (cckd->cdevhdr[sfx].nullfmt > CKDDASD_NULLTRK_FMTMAX)
        cckd->cdevhdr[sfx].nullfmt = 0;

    if (cckd->cdevhdr[sfx].nullfmt == 0 && dev->oslinux && dev->devtype == 0x3390)
        cckd->cdevhdr[sfx].nullfmt = CKDDASD_NULLTRK_FMT2;

    if (cckd->cdevhdr[sfx].nullfmt == CKDDASD_NULLTRK_FMT2)
        dev->oslinux = 1;

    return 0;

} /* end function cckd_read_chdr */

/*-------------------------------------------------------------------*/
/* Write compressed dasd header                                      */
/*-------------------------------------------------------------------*/
int cckd_write_chdr (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* File index                */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "file[%d] write_chdr\n", sfx);

    /* Set version.release.modlvl */
    cckd->cdevhdr[sfx].vrm[0] = CCKD_VERSION;
    cckd->cdevhdr[sfx].vrm[1] = CCKD_RELEASE;
    cckd->cdevhdr[sfx].vrm[2] = CCKD_MODLVL;

    if (cckd_write (dev, sfx, CCKD_DEVHDR_POS, &cckd->cdevhdr[sfx], CCKD_DEVHDR_SIZE) < 0)
        return -1;

    return 0;

} /* end function cckd_write_chdr */

/*-------------------------------------------------------------------*/
/* Read the level 1 table                                            */
/*-------------------------------------------------------------------*/
int cckd_read_l1 (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* File index                */
int             len;                    /* Length of level 1 table   */
int             i;                      /* Work integer              */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "file[%d] read_l1 offset 0x%"I64_FMT"x\n",
                sfx, (U64)CCKD_L1TAB_POS);

    /* Free the old level 1 table if it exists */
    cckd->l1[sfx] = cckd_free (dev, "l1", cckd->l1[sfx]);

    /* Allocate the level 1 table */
    len = cckd->cdevhdr[sfx].numl1tab * CCKD_L1ENT_SIZE;
    if ((cckd->l1[sfx] = cckd_malloc (dev, "l1", len)) == NULL)
        return -1;
    memset(cckd->l1[sfx], sfx ? 0xFF : 0, len);

    /* Read the level 1 table */
    if (cckd_read (dev, sfx, CCKD_L1TAB_POS, cckd->l1[sfx], len) < 0)
        return -1;

    /* Fix endianess */
    if (cckd->swapend[sfx])
        cckd_swapend_l1 (cckd->l1[sfx], cckd->cdevhdr[sfx].numl1tab);

    /* Determine bounds */
    cckd->l2bounds = CCKD_L1TAB_POS + len;
    for (i = 0; i < cckd->cdevhdr[sfx].numl1tab; i++)
        if (cckd->l1[sfx][i] != 0 && cckd->l1[sfx][i] != 0xffffffff)
            cckd->l2bounds += CCKD_L2TAB_SIZE;

    /* Check if all l2 tables are within bounds */
    cckd->l2ok = 1;
    for (i = 0; i < cckd->cdevhdr[sfx].numl1tab && cckd->l2ok; i++)
        if (cckd->l1[sfx][i] != 0 && cckd->l1[sfx][i] != 0xffffffff)
            if (cckd->l1[sfx][i] > cckd->l2bounds - CCKD_L2TAB_SIZE)
                cckd->l2ok = 0;

    return 0;

} /* end function cckd_read_l1 */

/*-------------------------------------------------------------------*/
/* Write the level 1 table                                           */
/*-------------------------------------------------------------------*/
int cckd_write_l1 (DEVBLK *dev)
{
CCKDDASD_EXT    *cckd;                  /* -> cckd extension         */
int             sfx;                    /* File index                */
int             len;                    /* Length of level 1 table   */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;
    len = cckd->cdevhdr[sfx].numl1tab * CCKD_L1ENT_SIZE;

    cckd_trace (dev, "file[%d] write_l1 0x%"I64_FMT"x len %d\n",
                sfx, (U64)CCKD_L1TAB_POS, len);

    if (cckd_write (dev, sfx, CCKD_L1TAB_POS, cckd->l1[sfx], len) < 0)
        return -1;

    return 0;

} /* end function cckd_write_l1 */

/*-------------------------------------------------------------------*/
/* Update a level 1 table entry                                      */
/*-------------------------------------------------------------------*/
int cckd_write_l1ent (DEVBLK *dev, int l1x)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* File index                */
off_t           off;                    /* Offset to l1 entry        */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;
    off = (off_t)(CCKD_L1TAB_POS + l1x * CCKD_L1ENT_SIZE);

    cckd_trace (dev, "file[%d] write_l1ent[%d] , 0x%" I64_FMT "x\n",
                sfx, l1x, (long long)off);

    if (cckd_write (dev, sfx, off, &cckd->l1[sfx][l1x], CCKD_L1ENT_SIZE) < 0)
        return -1;

    return 0;

} /* end function cckd_write_l1ent */

/*-------------------------------------------------------------------*/
/* Initial read                                                      */
/*-------------------------------------------------------------------*/
int cckd_read_init (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* File index                */
CKDDASD_DEVHDR  devhdr;                 /* Device header             */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "file[%d] read_init\n", sfx);

    /* Read the device header */
    if (cckd_read (dev, sfx, 0, &devhdr, CKDDASD_DEVHDR_SIZE) < 0)
        return -1;

    /* Check the device hdr */
    if (sfx == 0 && memcmp (&devhdr.devid, "CKD_C370", 8) == 0)
        cckd->ckddasd = 1;
    else if (sfx == 0 && memcmp (&devhdr.devid, "FBA_C370", 8) == 0)
        cckd->fbadasd = 1;
    else if (!(sfx && memcmp (&devhdr.devid, "CKD_S370", 8) == 0 && cckd->ckddasd)
          && !(sfx && memcmp (&devhdr.devid, "FBA_S370", 8) == 0 && cckd->fbadasd))
    {
        logmsg (_("HHCCD110E %4.4X file[%d] devhdr id error\n"),
                dev->devnum, sfx);
        return -1;
    }

    /* Read the compressed header */
    if (cckd_read_chdr (dev) < 0)
        return -1;

    /* Read the level 1 table */
    if (cckd_read_l1 (dev) < 0)
        return -1;

    return 0;
} /* end function cckd_read_init */

/*-------------------------------------------------------------------*/
/* Read free space                                                   */
/*-------------------------------------------------------------------*/
int cckd_read_fsp (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
off_t           fpos;                   /* Free space offset         */
int             sfx;                    /* File index                */
int             i;                      /* Index                     */
CCKD_FREEBLK    freeblk;                /* First freeblk read        */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    cckd_trace (dev, "file[%d] read_fsp number %d\n",
                sfx, cckd->cdevhdr[sfx].free_number);

    cckd->free = cckd_free (dev, "free", cckd->free);
    cckd->free1st = cckd->freelast = cckd->freeavail = -1;

    /* Get storage for the internal free space chain
     * in a multiple of 1024 entries
     */
    cckd->freenbr = (cckd->cdevhdr[sfx].free_number + 1023) & ~0x3FF;
    if (cckd->freenbr)
        if ((cckd->free = cckd_calloc (dev, "free", cckd->freenbr, CCKD_FREEBLK_ISIZE)) == NULL)
            return -1;

    /* Build the doubly linked internal free space chain */
    if (cckd->cdevhdr[sfx].free_number)
    {
        cckd->free1st = 0;

        /* Read the first freeblk to determine old/new format */
        fpos = (off_t)cckd->cdevhdr[sfx].free;
        if (cckd_read (dev, sfx, fpos, &freeblk, CCKD_FREEBLK_SIZE) < 0)
            return -1;

        if (memcmp(&freeblk, "FREE_BLK", 8) == 0)
        {
            /* new format free space */
            CCKD_FREEBLK *fsp;
            U32 ofree = cckd->cdevhdr[sfx].free;
            int n = cckd->cdevhdr[sfx].free_number * CCKD_FREEBLK_SIZE;
            if ((fsp = cckd_malloc (dev, "fsp", n)) == NULL)
                return -1;
            fpos += CCKD_FREEBLK_SIZE;
            if (cckd_read (dev, sfx, fpos, fsp, n) < 0)
                return -1;
            for (i = 0; i < cckd->cdevhdr[sfx].free_number; i++)
            {
                if (i == 0)
                    cckd->cdevhdr[sfx].free = fsp[i].pos;
                else
                    cckd->free[i-1].pos = fsp[i].pos;
                cckd->free[i].pos  = 0;
                cckd->free[i].len  = fsp[i].len;
                cckd->free[i].prev = i - 1;
                cckd->free[i].next = i + 1;
            }
            cckd->free[i-1].next = -1;
            cckd->freelast = i-1;
            fsp = cckd_free (dev, "fsp", fsp);

            /* truncate if new format free space was at the end */
            if (ofree == cckd->cdevhdr[sfx].size)
            {
                fpos = (off_t)cckd->cdevhdr[sfx].size;
                cckd_ftruncate(dev, sfx, fpos);
            }
        } /* new format free space */
        else
        {
            /* old format free space */
            fpos = (off_t)cckd->cdevhdr[sfx].free;
            for (i = 0; i < cckd->cdevhdr[sfx].free_number; i++)
            {
                if (cckd_read (dev, sfx, fpos, &cckd->free[i], CCKD_FREEBLK_SIZE) < 0)
                    return -1;
                cckd->free[i].prev = i - 1;
                cckd->free[i].next = i + 1;
                fpos = (off_t)cckd->free[i].pos;
            }
            cckd->free[i-1].next = -1;
            cckd->freelast = i-1;
        } /* old format free space */
    } /* if (cckd->cdevhdr[sfx].free_number) */

    /* Build singly linked chain of available free space entries */
    if (cckd->cdevhdr[sfx].free_number < cckd->freenbr)
    {
        cckd->freeavail = cckd->cdevhdr[sfx].free_number;
        for (i = cckd->freeavail; i < cckd->freenbr; i++)
            cckd->free[i].next = i + 1;
        cckd->free[i-1].next = -1;
    }

    /* Set minimum free space size */
    cckd->freemin = CCKD_FREE_MIN_SIZE + ((cckd->freenbr >> 10) * CCKD_FREE_MIN_INCR);
    return 0;

} /* end function cckd_read_fsp */

/*-------------------------------------------------------------------*/
/* Write the free space                                              */
/*-------------------------------------------------------------------*/
int cckd_write_fsp (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
off_t           fpos;                   /* Free space offset         */
U32             ppos;                   /* Previous free space offset*/
int             sfx;                    /* File index                */
int             i, j, n;                /* Work variables            */
int             rc;                     /* Return code               */
CCKD_FREEBLK   *fsp = NULL;             /* -> new format free space  */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;

    if (cckd->free == NULL)
        return 0;

    cckd_trace (dev, "file[%d] write_fsp number %d\n",
                sfx, cckd->cdevhdr[sfx].free_number);

    /* get rid of pending free space */
    for (i = 0; i < CCKD_MAX_FREEPEND; i++)
        cckd_flush_space(dev);

    /* sanity checks */
    if (cckd->cdevhdr[sfx].free_number == 0 || cckd->cdevhdr[sfx].free == 0)
    {
        cckd->cdevhdr[sfx].free_number = cckd->cdevhdr[sfx].free = 0;
        cckd->free1st = cckd->freelast = cckd->freeavail = -1;
    }

    /* Write any free spaces */
    if (cckd->cdevhdr[sfx].free)
    {
        /* size needed for new format free space */
        n = (cckd->cdevhdr[sfx].free_number+1) * CCKD_FREEBLK_SIZE;

        /* look for existing free space to fit new format free space */   
        fpos = 0;
        for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
            if (n <= (int)cckd->free[i].len)
                break;
        if (i >= 0)
            fpos = cckd->free[i].prev < 0
                 ? (off_t)cckd->cdevhdr[sfx].free
                 : (off_t)cckd->free[cckd->free[i].prev].pos;

        /* if no applicable space see if we can append to the file */
        if (fpos == 0 && cckd->maxsize - cckd->cdevhdr[sfx].size >= n)
            fpos = (off_t)cckd->cdevhdr[sfx].size;

        if (fpos && (fsp = cckd_malloc (dev, "fsp", n)) == NULL)
            fpos = 0;

        if (fpos)
        {
            /* New format free space */
            memcpy (&fsp[0], "FREE_BLK", 8);
            ppos = cckd->cdevhdr[sfx].free;
            for (i = cckd->free1st, j = 1; i >= 0; i = cckd->free[i].next)
            {
                fsp[j].pos = ppos;
                fsp[j++].len = cckd->free[i].len;
                ppos = cckd->free[i].pos;
            }
            rc = cckd_write (dev, sfx, fpos, fsp, n);
            fsp = cckd_free (dev, "fsp", fsp);
            if (rc < 0)
                return -1;
            cckd->cdevhdr[sfx].free = (U32)fpos;
        } /* new format free space */
        else
        {
            /* Old format free space */
            for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
            {
                if (cckd_write (dev, sfx, fpos, &cckd->free[i], CCKD_FREEBLK_SIZE) < 0)
                    return -1;
                fpos = (off_t)cckd->free[i].pos;
            }
        } /* old format free space */
    } /* if (cckd->cdevhdr[sfx].free) */

    /* Free the free space array */
    cckd->free = cckd_free (dev, "free", cckd->free);
    cckd->freenbr = 0;
    cckd->free1st = cckd->freelast = cckd->freeavail = -1;

    return 0;

} /* end function cckd_write_fsp */

/*-------------------------------------------------------------------*/
/* Read a new level 2 table                                          */
/*-------------------------------------------------------------------*/
int cckd_read_l2 (DEVBLK *dev, int sfx, int l1x)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
off_t           off;                    /* L2 file offset            */
int             fnd;                    /* Found cache               */
int             lru;                    /* Oldest available cache    */
CCKD_L2ENT     *buf;                    /* -> Cache buffer           */
int             i;                      /* Loop index                */
int             nullfmt;                /* Null track format         */

    cckd = dev->cckd_ext;
    nullfmt = cckd->cdevhdr[cckd->sfn].nullfmt;

    cckd_trace (dev, "file[%d] read_l2 %d active %d %d %d\n",
                sfx, l1x, cckd->sfx, cckd->l1x, cckd->l2active);

    /* Return if table is already active */
    if (sfx == cckd->sfx && l1x == cckd->l1x) return 0;

    cache_lock(CACHE_L2);

    /* Inactivate the previous entry */
    if (cckd->l2active >= 0)
        cache_setflag(CACHE_L2, cckd->l2active, ~L2_CACHE_ACTIVE, 0);
    cckd->l2 = NULL;
    cckd->l2active = cckd->sfx = cckd->l1x = -1;

    /* scan the cache array for the l2tab */
    fnd = cache_lookup (CACHE_L2, L2_CACHE_SETKEY(sfx, dev->devnum, l1x), &lru);

    /* check for level 2 cache hit */
    if (fnd >= 0)
    {
        cckd_trace (dev, "l2[%d,%d] cache[%d] hit\n", sfx, l1x, fnd);
        cache_setflag (CACHE_L2, fnd, 0, L2_CACHE_ACTIVE);
        cache_setage (CACHE_L2, fnd);
        cckdblk.stats_l2cachehits++;
        cache_unlock (CACHE_L2);
        cckd->sfx = sfx;
        cckd->l1x = l1x;
        cckd->l2 = cache_getbuf(CACHE_L2, fnd, 0);
        cckd->l2active = fnd;
        return 1;
    }

    cckd_trace (dev, "l2[%d,%d] cache[%d] miss\n", sfx, l1x, lru);

    /* Steal an entry if all are busy */
    if (lru < 0) lru = cckd_steal_l2();

    /* Make the entry active */
    cache_setkey (CACHE_L2, lru, L2_CACHE_SETKEY(sfx, dev->devnum, l1x));
    cache_setflag (CACHE_L2, lru, 0, L2_CACHE_ACTIVE);
    cache_setage (CACHE_L2, lru);
    buf = cache_getbuf(CACHE_L2, lru, CCKD_L2TAB_SIZE);
    cckdblk.stats_l2cachemisses++;
    cache_unlock (CACHE_L2);
    if (buf == NULL) return -1;

    /* Check for null table */
    if (cckd->l1[sfx][l1x] == 0)
    {
        memset(buf, 0, CCKD_L2TAB_SIZE);
        if (nullfmt)
            for (i = 0; i < 256; i++)
                buf[i].len = buf[i].size = nullfmt;
        cckd_trace (dev, "l2[%d,%d] cache[%d] null fmt[%d]\n", sfx, l1x, lru, nullfmt);
    }
    else if (cckd->l1[sfx][l1x] == 0xffffffff)
    {
        memset(buf, 0xff, CCKD_L2TAB_SIZE);
        cckd_trace (dev, "l2[%d,%d] cache[%d] null 0xff\n", sfx, l1x, lru);
    }
    /* Read the new level 2 table */
    else
    {
        off = (off_t)cckd->l1[sfx][l1x];
        if (cckd_read (dev, sfx, off, buf, CCKD_L2TAB_SIZE) < 0)
        {
            cache_lock(CACHE_L2);
            cache_setflag(CACHE_L2, lru, 0, 0);
            cache_unlock(CACHE_L2);
            return -1;
        }

        if (cckd->swapend[sfx])
            cckd_swapend_l2 (buf);

        cckd_trace (dev, "file[%d] cache[%d] l2[%d] read offset 0x%" I64_FMT "x\n",
                    sfx, lru, l1x, (long long)cckd->l1[sfx][l1x]);

        cckd->l2reads[sfx]++;
        cckd->totl2reads++;
        cckdblk.stats_l2reads++;
    }

    cckd->sfx = sfx;
    cckd->l1x = l1x;
    cckd->l2 = buf;
    cckd->l2active = lru;

    return 0;

} /* end function cckd_read_l2 */

/*-------------------------------------------------------------------*/
/* Purge all l2tab cache entries for a given device                  */
/*-------------------------------------------------------------------*/
void cckd_purge_l2 (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "purge_l2%s\n", "");

    cache_lock (CACHE_L2);
    cckd->l2active = cckd->sfx = cckd->l1x = -1;
    cckd->l2 = NULL;
    cache_scan (CACHE_L2, cckd_purge_l2_scan, dev);
    cache_unlock (CACHE_L2);
}
int cckd_purge_l2_scan (int *answer, int ix, int i, void *data)
{
U16             sfx;                    /* Cached suffix             */
U16             devnum;                 /* Cached device number      */
U32             l1x;                    /* Cached level 1 index      */
DEVBLK         *dev = data;             /* -> device block           */

    UNREFERENCED(answer);
    L2_CACHE_GETKEY(i, sfx, devnum, l1x);
    if (dev == NULL || devnum == dev->devnum)
    {
        cckd_trace (dev, "purge l2cache[%d] %4.4X sfx %d ix %d purged\n",
                    i, devnum, sfx, l1x);
        cache_release(ix, i, 0);
    }
    return 0;
}

/*-------------------------------------------------------------------*/
/* Steal an l2tab cache entry                                        */
/*-------------------------------------------------------------------*/
int cckd_steal_l2 ()
{
DEVBLK         *dev;                    /* -> device block           */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             i;                      /* Stolen cache index        */
U16             sfx;                    /* Cached suffix             */
U16             devnum;                 /* Cached device number      */
U32             l1x;                    /* Cached level 1 index      */

    i = cache_scan (CACHE_L2, cckd_steal_l2_scan, NULL);
    L2_CACHE_GETKEY(i, sfx, devnum, l1x);
    dev = cckd_find_device_by_devnum(devnum);
    cckd = dev->cckd_ext;
    cckd->l2active = cckd->sfx = cckd->l1x = -1;
    cckd->l2 = NULL;
    cache_release(CACHE_L2, i, 0);
    return i;
}
int cckd_steal_l2_scan (int *answer, int ix, int i, void *data)
{
    UNREFERENCED(data);
    if (*answer < 0) *answer = i;
    else if (cache_getage(ix, i) < cache_getage(ix, *answer))
        *answer = i;
    return 0;
}

/*-------------------------------------------------------------------*/
/* Write the current level 2 table                                   */
/*-------------------------------------------------------------------*/
int cckd_write_l2 (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx,l1x;                /* Lookup table indices      */
off_t           off, old_off;           /* New/old L2 file offsets   */
int             size = CCKD_L2TAB_SIZE; /* L2 table size             */
int             fix;                    /* Null format type          */

    cckd = dev->cckd_ext;
    sfx = cckd->sfn;
    l1x = cckd->l1x;
    fix = cckd->cdevhdr[sfx].nullfmt;
    cckd->l2ok = 0;

    cckd_trace (dev, "file[%d] write_l2 %d\n", sfx, l1x);

    if (sfx < 0 || l1x < 0) return -1;

    old_off = (off_t)cckd->l1[sfx][l1x];

    if (cckd->l1[sfx][l1x] == 0 || cckd->l1[sfx][l1x] == 0xffffffff)
        cckd->l2bounds += CCKD_L2TAB_SIZE;

    /* Write the L2 table if it's not empty */
    if (memcmp(cckd->l2, &empty_l2[fix], CCKD_L2TAB_SIZE))
    {
        if ((off = cckd_get_space (dev, &size, CCKD_L2SPACE)) < 0)
            return -1;
        if (cckd_write (dev, sfx, off, cckd->l2, CCKD_L2TAB_SIZE) < 0)
            return -1;
    }
    else
    {
        off = 0;
        cckd->l2bounds -= CCKD_L2TAB_SIZE;
    }

    /* Free the old L2 space */
    cckd_rel_space (dev, old_off, CCKD_L2TAB_SIZE, CCKD_L2TAB_SIZE);

    /* Update level 1 table */
    cckd->l1[sfx][l1x] = (U32)off;
    if (cckd_write_l1ent (dev, l1x) < 0)
        return -1;

    return 0;

} /* end function cckd_write_l2 */

/*-------------------------------------------------------------------*/
/* Return a level 2 entry                                            */
/*-------------------------------------------------------------------*/
int cckd_read_l2ent (DEVBLK *dev, CCKD_L2ENT *l2, int trk)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx,l1x,l2x;            /* Lookup table indices      */

    cckd = dev->cckd_ext;

    l1x = trk >> 8;
    l2x = trk & 0xff;

    if (l2 != NULL) l2->pos = l2->len = l2->size = 0;

    for (sfx = cckd->sfn; sfx >= 0; sfx--)
    {
        cckd_trace (dev, "file[%d] l2[%d,%d] trk[%d] read_l2ent 0x%x\n",
                    sfx, l1x, l2x, trk, cckd->l1[sfx][l1x]);

        /* Continue if l2 table not in this file */
        if (cckd->l1[sfx][l1x] == 0xffffffff)
            continue;

        /* Read l2 table from this file */
        if (cckd_read_l2 (dev, sfx, l1x) < 0)
            return -1;

        /* Exit loop if track is in this file */
        if (cckd->l2[l2x].pos != 0xffffffff)
            break;
    }

    cckd_trace (dev, "file[%d] l2[%d,%d] trk[%d] read_l2ent 0x%x %d %d\n",
                sfx, l1x, l2x, trk, sfx >= 0 ? cckd->l2[l2x].pos : 0,
                sfx >= 0 ? cckd->l2[l2x].len : 0,
                sfx >= 0 ? cckd->l2[l2x].size : 0);

    if (l2 != NULL && sfx >= 0)
    {
        l2->pos = cckd->l2[l2x].pos;
        l2->len = cckd->l2[l2x].len;
        l2->size = cckd->l2[l2x].size;
    }

    return sfx;

} /* end function cckd_read_l2ent */

/*-------------------------------------------------------------------*/
/* Update a level 2 entry                                            */
/*-------------------------------------------------------------------*/
int cckd_write_l2ent (DEVBLK *dev,  CCKD_L2ENT *l2, int trk)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx,l1x,l2x;            /* Lookup table indices      */
off_t           off;                    /* L2 entry offset           */

    cckd = dev->cckd_ext;

    /* Error return if no available level 2 table */
    if (!cckd->l2) return -1;

    sfx = cckd->sfn;
    l1x = trk >> 8;
    l2x = trk & 0xff;

    /* Copy the new entry if passed */
    if (l2) memcpy (&cckd->l2[l2x], l2, CCKD_L2ENT_SIZE);

    cckd_trace (dev, "file[%d] l2[%d,%d] trk[%d] write_l2ent 0x%x %d %d\n",
                sfx, l1x, l2x, trk,
                cckd->l2[l2x].pos, cckd->l2[l2x].len, cckd->l2[l2x].size);

    /* If no level 2 table for this file, then write a new one */
    if (cckd->l1[sfx][l1x] == 0 || cckd->l1[sfx][l1x] == 0xffffffff)
        return cckd_write_l2 (dev);

    /* Write the level 2 table entry */
    off = (off_t)(cckd->l1[sfx][l1x] + l2x * CCKD_L2ENT_SIZE);
    if (cckd_write (dev, sfx, off, &cckd->l2[l2x], CCKD_L2ENT_SIZE) < 0)
        return -1;

    return 0;
} /* end function cckd_write_l2ent */

/*-------------------------------------------------------------------*/
/* Read a track image                                                */
/*-------------------------------------------------------------------*/
int cckd_read_trkimg (DEVBLK *dev, BYTE *buf, int trk, BYTE *unitstat)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             sfx;                    /* File index                */
CCKD_L2ENT      l2;                     /* Level 2 entry             */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "trk[%d] read_trkimg\n", trk);

    /* Read level 2 entry for the track */
    if ((sfx = cckd_read_l2ent (dev, &l2, trk)) < 0)
        goto cckd_read_trkimg_error;

    /* Read the track image or build a null track image */
    if (l2.pos != 0)
    {
        rc = cckd_read (dev, sfx, (off_t)l2.pos, buf, (size_t)l2.len);
        if (rc < 0)
            goto cckd_read_trkimg_error;

        cckd->reads[sfx]++;
        cckd->totreads++;
        cckdblk.stats_reads++;
        cckdblk.stats_readbytes += rc;
        if (cckd->notnull == 0 && trk > 1) cckd->notnull = 1;
    }
    else
        rc = cckd_null_trk (dev, buf, trk, l2.len);

    /* Validate the track image */
    if (cckd_cchh (dev, buf, trk) < 0)
        goto cckd_read_trkimg_error;

    return rc;

cckd_read_trkimg_error:

    if (unitstat)
    {
        ckd_build_sense (dev, SENSE_EC, 0, 0, FORMAT_1, MESSAGE_0);
        *unitstat = CSW_CE | CSW_DE | CSW_UC;
    }

    return cckd_null_trk (dev, buf, trk, 0);

} /* end function cckd_read_trkimg */

/*-------------------------------------------------------------------*/
/* Write a track image                                               */
/*-------------------------------------------------------------------*/
int cckd_write_trkimg (DEVBLK *dev, BYTE *buf, int len, int trk, int flags)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
off_t           off;                    /* File offset               */
CCKD_L2ENT      l2, oldl2;              /* Level 2 entries           */
int             sfx,l1x,l2x;            /* Lookup table indices      */
int             after = 0;              /* 1=New track after old     */
int             size;                   /* Size of new track         */

    cckd = dev->cckd_ext;

    sfx = cckd->sfn;
    l1x = trk >> 8;
    l2x = trk & 0xff;

    cckd_trace (dev, "file[%d] trk[%d] write_trkimg len %d buf %p:%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                sfx, trk, len, buf, buf[0], buf[1], buf[2], buf[3], buf[4]);

    /* Validate the new track image */
    if (cckd_cchh (dev, buf, trk) < 0)
        return -1;

    /* Get the level 2 table for the track in the active file */
    if (cckd_read_l2 (dev, sfx, l1x) < 0)
        return -1;

    /* Save the level 2 entry for the track */
    oldl2.pos = cckd->l2[l2x].pos;
    oldl2.len = cckd->l2[l2x].len;
    oldl2.size = cckd->l2[l2x].size;
    cckd_trace (dev, "file[%d] trk[%d] write_trkimg oldl2 0x%x %d %d\n",
                sfx, trk, oldl2.pos,oldl2.len,oldl2.size);

    /* Check if writing a null track */
    len = cckd_check_null_trk(dev, buf, trk, len);

    if (len > CKDDASD_NULLTRK_FMTMAX)
    {
        /* Get space for the track image */
        size = len;
        if ((off = cckd_get_space (dev, &size, flags)) < 0)
            return -1;

        l2.pos = (U32)off;
        l2.len = (U16)len;
        l2.size = (U16)size;

        if (oldl2.pos != 0 && oldl2.pos != 0xffffffff && oldl2.pos < l2.pos)
            after = 1;

        /* Write the track image */
        if ((rc = cckd_write (dev, sfx, off, buf, len)) < 0)
            return -1;

        cckd->writes[sfx]++;
        cckd->totwrites++;
        cckdblk.stats_writes++;
        cckdblk.stats_writebytes += rc;
    }
    else
    {
        l2.pos = 0;
        l2.len = l2.size = (U16)len;
    }

    /* Update the level 2 entry */
    if (cckd_write_l2ent (dev, &l2, trk) < 0)
        return -1;

    /* Release the previous space */
    cckd_rel_space (dev, (off_t)oldl2.pos, (int)oldl2.len, (int)oldl2.size);

    /* `after' is 1 if the new offset is after the old offset */
    return after;

} /* end function cckd_write_trkimg */

/*-------------------------------------------------------------------*/
/* Harden the file                                                   */
/*-------------------------------------------------------------------*/
int cckd_harden(DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc=0;                   /* Return code               */

    cckd = dev->cckd_ext;

    if ((dev->ckdrdonly && cckd->sfn == 0)
     || cckd->open[cckd->sfn] != CCKD_OPEN_RW)
        return 0;

    cckd_trace (dev, "file[%d] harden\n", cckd->sfn);

    /* Write the compressed device header */
    if (cckd_write_chdr (dev) < 0)
        rc = -1;

    /* Write the level 1 table */
    if (cckd_write_l1 (dev) < 0)
        rc = -1;

    /* Write the free space chain */
    if (cckd_write_fsp (dev) < 0)
        rc = -1;

    /* Re-write the compressed device header */
    cckd->cdevhdr[cckd->sfn].options &= ~CCKD_OPENED;
    if (cckd_write_chdr (dev) < 0)
        rc = -1;

    if (cckdblk.fsync)
        fdatasync (cckd->fd[cckd->sfn]);

    return rc;
} /* cckd_harden */

/*-------------------------------------------------------------------*/
/* Return length of an uncompressed track image                      */
/*-------------------------------------------------------------------*/
int cckd_trklen (DEVBLK *dev, BYTE *buf)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             size;                   /* Track size                */

    cckd = dev->cckd_ext;

    if (cckd->fbadasd)
        return CKDDASD_TRKHDR_SIZE + CFBA_BLOCK_SIZE;

    for (size = CKDDASD_TRKHDR_SIZE;
         memcmp (&buf[size], &eighthexFF, 8) != 0; )
    {
        if (size > dev->ckdtrksz) break;

        /* add length of count, key, and data fields */
        size += CKDDASD_RECHDR_SIZE +
                buf[size+5] +
                (buf[size+6] << 8) + buf[size+7];
    }

    /* add length for end-of-track indicator */
    size += CKDDASD_RECHDR_SIZE;

    /* check for missing end-of-track indicator */
    if (size > dev->ckdtrksz ||
        memcmp (&buf[size-CKDDASD_RECHDR_SIZE], &eighthexFF, 8) != 0)
    {
        logmsg (_("HHCCD121E %4.4X file[%d] trklen err for %2.2x%2.2x%2.2x%2.2x%2.2x\n"),
                dev->devnum, cckd->sfn, buf[0], buf[1], buf[2], buf[3], buf[4]);
        size = -1;
    }

    return size;
}

/*-------------------------------------------------------------------*/
/* Build a null track                                                */
/*-------------------------------------------------------------------*/
int cckd_null_trk(DEVBLK *dev, BYTE *buf, int trk, int nullfmt)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             i;                      /* Loop counter              */
CKDDASD_TRKHDR *trkhdr;                 /* -> Track header           */
CKDDASD_RECHDR *rechdr;                 /* -> Record header          */
U32             cyl;                    /* Cylinder number           */
U32             head;                   /* Head number               */
BYTE            r;                      /* Record number             */
BYTE           *pos;                    /* -> Next position in buffer*/
int             len;                    /* Length of null track      */

    cckd = dev->cckd_ext;

    if (nullfmt < 0 || nullfmt > CKDDASD_NULLTRK_FMTMAX)
        nullfmt = cckd->cdevhdr[cckd->sfn].nullfmt;

    // FIXME
    // Compatibility check for nullfmt bug and linux -- 18 May 2005
    // Remove at some reasonable date in the future
    else if (nullfmt == 0
     && cckd->cdevhdr[cckd->sfn].nullfmt == CKDDASD_NULLTRK_FMT2)
        nullfmt = CKDDASD_NULLTRK_FMT2;

    if (cckd->ckddasd)
    {

        /* cylinder and head calculations */
        cyl = trk / dev->ckdheads;
        head = trk % dev->ckdheads;

        /* Build the track header */
        trkhdr = (CKDDASD_TRKHDR*)buf;
        trkhdr->bin = 0;
        store_hw(&trkhdr->cyl, cyl);
        store_hw(&trkhdr->head, head);
        pos = buf + CKDDASD_TRKHDR_SIZE;

        /* Build record zero */
        r = 0;
        rechdr = (CKDDASD_RECHDR*)pos;
        pos += CKDDASD_RECHDR_SIZE;
        store_hw(&rechdr->cyl, cyl);
        store_hw(&rechdr->head, head);
        rechdr->rec = r;
        rechdr->klen = 0;
        store_hw(&rechdr->dlen, 8);
        memset (pos, 0, 8);
        pos += 8;
        r++;

        /* Specific null track formatting */
        if (nullfmt == CKDDASD_NULLTRK_FMT0)
        {
            rechdr = (CKDDASD_RECHDR*)pos;
            pos += CKDDASD_RECHDR_SIZE;

            store_hw(&rechdr->cyl, cyl);
            store_hw(&rechdr->head, head);
            rechdr->rec = r;
            rechdr->klen = 0;
            store_hw(&rechdr->dlen, 0);
            r++;
        }
        else if (nullfmt == CKDDASD_NULLTRK_FMT2)
        {
            for (i = 0; i < 12; i++)
            {
                rechdr = (CKDDASD_RECHDR*)pos;
                pos += CKDDASD_RECHDR_SIZE;

                store_hw(&rechdr->cyl, cyl);
                store_hw(&rechdr->head, head);
                rechdr->rec = r;
                rechdr->klen = 0;
                store_hw(&rechdr->dlen, 4096);
                r++;
                memset(pos, 0, 4096);
                pos += 4096;
            }
        }

        /* Build the end of track marker */
        memcpy (pos, eighthexFF, 8);
        pos += 8;
        len = (int)(pos - buf);
    }
    else
    {
        memset (buf, 0, CFBA_BLOCK_SIZE + CKDDASD_TRKHDR_SIZE);
        store_fw(buf+1, trk);
        len = CFBA_BLOCK_SIZE + CKDDASD_TRKHDR_SIZE;
    }

    cckd_trace (dev, "null_trk %s %d format %d size %d\n",
                cckd->ckddasd ? "trk" : "blkgrp", trk, nullfmt, len);

    return len;

} /* end function cckd_null_trk */

/*-------------------------------------------------------------------*/
/* Return a number 0 .. CKDDASD_NULLTRK_FMTMAX if track is null      */
/* else return the original length                                   */
/*-------------------------------------------------------------------*/
int cckd_check_null_trk (DEVBLK *dev, BYTE *buf, int trk, int len)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
BYTE            buf2[65536];            /* Null track buffer         */

    cckd = dev->cckd_ext;
    rc = len;

    if (len == CKDDASD_NULLTRK_SIZE0)
        rc = CKDDASD_NULLTRK_FMT0;
    else if (len == CKDDASD_NULLTRK_SIZE1)
        rc = CKDDASD_NULLTRK_FMT1;
    else if (len == CKDDASD_NULLTRK_SIZE2 && dev->oslinux
          && (!cckd->notnull || cckdblk.linuxnull))
    {
         cckd_null_trk (dev, buf2, trk, 0);
         if (memcmp(buf, buf2, len) == 0)
            rc = CKDDASD_NULLTRK_FMT2;
    }

    return rc;
}

/*-------------------------------------------------------------------*/
/* Verify a track/block header and return track/block number         */
/*-------------------------------------------------------------------*/
int cckd_cchh (DEVBLK *dev, BYTE *buf, int trk)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
U16             cyl;                    /* Cylinder                  */
U16             head;                   /* Head                      */
int             t;                      /* Calculated track          */
BYTE            badcomp=0;              /* 1=Unsupported compression */
static char    *comp[] = {"none", "zlib", "bzip2"};

    cckd = dev->cckd_ext;

    /* CKD dasd header verification */
    if (cckd->ckddasd)
    {
        cyl = fetch_hw (buf + 1);
        head = fetch_hw (buf + 3);
        t = cyl * dev->ckdheads + head;

        if (cyl < dev->ckdcyls && head < dev->ckdheads
         && (trk == -1 || t == trk))
        {
            if (buf[0] & ~cckdblk.comps)
            {
                if (buf[0] & ~CCKD_COMPRESS_MASK)
                {
                    if (cckdblk.bytemsgs++ < 10)
                        logmsg (_("HHCCD122E %4.4X file[%d] invalid byte 0 trk %d: "
                            "buf %2.2x%2.2x%2.2x%2.2x%2.2x\n"), dev->devnum, cckd->sfn,
                            t, buf[0],buf[1],buf[2],buf[3],buf[4]);
                    buf[0] &= CCKD_COMPRESS_MASK;
                }
            }
            if (buf[0] & ~cckdblk.comps)
                badcomp = 1;
            else
                return t;
        }
    }
    /* FBA dasd header verification */
    else
    {
        t = fetch_fw (buf + 1);
        if (t < dev->fbanumblk && (trk == -1 || t == trk))
        {
            if (buf[0] & ~cckdblk.comps)
            {
                if (buf[0] & ~CCKD_COMPRESS_MASK)
                {
                    logmsg (_("HHCCD123E %4.4X file[%d] invalid byte 0 blkgrp %d: "
                            "buf %2.2x%2.2x%2.2x%2.2x%2.2x\n"), dev->devnum, cckd->sfn,
                            t, buf[0],buf[1],buf[2],buf[3],buf[4]);
                    buf[0] &= CCKD_COMPRESS_MASK;
                }
            }
            if (buf[0] & ~cckdblk.comps)
                badcomp = 1;
            else
                return t;
        }
    }

    if (badcomp)
    {
        logmsg (_("HHCCD124E %4.4X file[%d] invalid %s hdr %s %d: "
                "%s compression unsupported\n"),
                dev->devnum, cckd->sfn,
                cckd->ckddasd ? "trk" : "blk",
                cckd->ckddasd ? "trk" : "blk", t, comp[buf[0]]);
    }
    else
    {
        logmsg (_("HHCCD125E %4.4X file[%d] invalid %s hdr %s %d "
                "buf %p:%2.2x%2.2x%2.2x%2.2x%2.2x\n"),
                dev->devnum, cckd->sfn,
                cckd->ckddasd ? "trk" : "blk",
                cckd->ckddasd ? "trk" : "blk", trk,
                buf, buf[0], buf[1], buf[2], buf[3], buf[4]);
        cckd_print_itrace ();
    }

    return -1;
} /* end function cckd_cchh */

/*-------------------------------------------------------------------*/
/* Validate a track image                                            */
/*-------------------------------------------------------------------*/
int cckd_validate (DEVBLK *dev, BYTE *buf, int trk, int len)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             cyl;                    /* Cylinder                  */
int             head;                   /* Head                      */
char            cchh[4],cchh2[4];       /* Cyl, head big-endian      */
int             r;                      /* Record number             */
int             sz;                     /* Track size                */
int             vlen;                   /* Validation length         */
int             kl,dl;                  /* Key/Data lengths          */

    cckd = dev->cckd_ext;

    if (buf == NULL || len < 0) return -1;

    cckd_trace (dev, "validating %s %d len %d %2.2x%2.2x%2.2x%2.2x%2.2x "
                "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                cckd->ckddasd ? "trk" : "blkgrp", trk, len,
                buf[0], buf[1], buf[2], buf[3], buf[4], buf[5], buf[6],
                buf[7], buf[8], buf[9], buf[10], buf[11], buf[12]);

    /* FBA dasd check */
    if (cckd->fbadasd)
    {
        if (len == CFBA_BLOCK_SIZE + CKDDASD_TRKHDR_SIZE || len == 0)
            return len;
        cckd_trace (dev, "validation failed: bad length%s\n","");
        return -1;
    }

    /* cylinder and head calculations */
    cyl = trk / dev->ckdheads;
    head = trk % dev->ckdheads;
    cchh[0] = cyl >> 8;
    cchh[1] = cyl & 0xFF;
    cchh[2] = head >> 8;
    cchh[3] = head & 0xFF;

    /* validate record 0 */
    memcpy (cchh2, &buf[5], 4); cchh2[0] &= 0x7f; /* fix for ovflow */
    if (/* memcmp (cchh, cchh2, 4) != 0 || */ buf[9]  != 0 ||
        buf[10] != 0 || buf[11] != 0 || buf[12] != 8)
    {
        cckd_trace (dev, "validation failed: bad r0%s\n","");
        return -1;
    }

    /* validate records 1 thru n */
    vlen = len > 0 ? len : dev->ckdtrksz;
    for (r = 1, sz = 21; sz + 8 <= vlen; sz += 8 + kl + dl, r++)
    {
        if (memcmp (&buf[sz], eighthexFF, 8) == 0) break;
        kl = buf[sz+5];
        dl = buf[sz+6] * 256 + buf[sz+7];
        /* fix for track overflow bit */
        memcpy (cchh2, &buf[sz], 4); cchh2[0] &= 0x7f;

        /* fix for funny formatted vm disks */
        /*
        if (r == 1) memcpy (cchh, cchh2, 4);
        */

        if (/*memcmp (cchh, cchh2, 4) != 0 ||*/ buf[sz+4] == 0 ||
            sz + 8 + kl + dl >= vlen)
        {
            cckd_trace (dev, "validation failed: bad r%d "
                        "%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x%2.2x\n",
                        r, buf[sz], buf[sz+1], buf[sz+2], buf[sz+3],
                        buf[sz+4], buf[sz+5], buf[sz+6], buf[sz+7]);
             return -1;
        }
    }
    sz += 8;

    if ((sz != len && len > 0) || sz > vlen)
    {
        cckd_trace (dev, "validation failed: no eot%s\n","");
        return -1;
    }

    return sz;

} /* end function cckd_validate */

/*-------------------------------------------------------------------*/
/* Return shadow file name                                           */
/*-------------------------------------------------------------------*/
char *cckd_sf_name (DEVBLK *dev, int sfx)
{
    /* Return base file name if index is 0 */
    if (sfx == 0)
        return dev->filename;

    /* Error if no shadow file name specified or number exceeded */
    if (dev->dasdsfn == NULL || sfx > CCKD_MAX_SF)
        return NULL;

    /* Set the suffix character in the shadow file name */
    if (sfx > 0)
        *dev->dasdsfx = '0' + sfx;
    else
        *dev->dasdsfx = '*';

    return dev->dasdsfn;

} /* end function cckd_sf_name */

/*-------------------------------------------------------------------*/
/* Initialize shadow files                                           */
/*-------------------------------------------------------------------*/
int cckd_sf_init (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
int             i;                      /* Index                     */
struct stat     st;                     /* stat() buffer             */
char            pathname[MAX_PATH];     /* file path in host format  */

    cckd = dev->cckd_ext;

    /* return if no shadow files */
    if (dev->dasdsfn == NULL) return 0;

#if 1
    /* Check for shadow file name collision */
    for (i = 1; i <= CCKD_MAX_SF && dev->dasdsfn != NULL; i++)
    {
     DEVBLK       *dev2;
     CCKDDASD_EXT *cckd2;
     int           j;

        for (dev2 = cckdblk.dev1st; dev2; dev2 = cckd2->devnext)
        {
            cckd2 = dev2->cckd_ext;
            if (dev2 == dev) continue;
            for (j = 0; j <= CCKD_MAX_SF && dev2->dasdsfn != NULL; j++)
            {
                if (strcmp (cckd_sf_name(dev, i),cckd_sf_name(dev2, j)) == 0)
                {
                    logmsg (_("HHCCD142E %4.4X file[%d] shadow file name %s\n"
                            "      collides with %4.4X file[%d] name %s\n"),
                            dev->devnum, i, cckd_sf_name(dev, i),
                            dev2->devnum, j, cckd_sf_name(dev2, j));
                    return -1;
                }
            }
        }
    }
#endif

    /* open all existing shadow files */
    for (cckd->sfn = 1; cckd->sfn <= CCKD_MAX_SF; cckd->sfn++)
    {
        hostpath(pathname, cckd_sf_name (dev, cckd->sfn), sizeof(pathname));
        if (stat (pathname, &st) < 0)
            break;

        /* Try to open the shadow file read-write then read-only */
        if (cckd_open (dev, cckd->sfn, O_RDWR|O_BINARY, 1) < 0)
            if (cckd_open (dev, cckd->sfn, O_RDONLY|O_BINARY, 0) < 0)
                break;

        /* Call the chkdsk function */
        rc = cckd_chkdsk (dev, 0);
        if (rc < 0) return -1;

        /* Perform initial read */
        rc = cckd_read_init (dev);
    }

    /* Backup to the last opened file number */
    cckd->sfn--;

    /* If the last file was opened read-only then create a new one   */
    if (cckd->open[cckd->sfn] == CCKD_OPEN_RO)
        if (cckd_sf_new(dev) < 0)
            return -1;

    /* Re-open previous rdwr files rdonly */
    for (i = 0; i < cckd->sfn; i++)
    {
        if (cckd->open[i] == CCKD_OPEN_RO) continue;
        if (cckd_open (dev, i, O_RDONLY|O_BINARY, 0) < 0)
        {
            logmsg (_("HHCCD151E %4.4X file[%d] error re-opening %s readonly\n  %s\n"),
                    dev->devnum, i, cckd_sf_name(dev, i), strerror(errno));
            return -1;
        }
    }

    return 0;

} /* end function cckd_sf_init */

/*-------------------------------------------------------------------*/
/* Create a new shadow file                                          */
/*-------------------------------------------------------------------*/
int cckd_sf_new (DEVBLK *dev)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             l1size;                 /* Size of level 1 table     */
CKDDASD_DEVHDR  devhdr;                 /* Device header             */

    cckd = dev->cckd_ext;

    cckd_trace (dev, "file[%d] sf_new %s\n", cckd->sfn+1,
                cckd_sf_name(dev, cckd->sfn+1) ?
                (char *)cckd_sf_name(dev, cckd->sfn+1) : "(none)");

    /* Error if no shadow file name */
    if (dev->dasdsfn == NULL)
    {
        logmsg (_("HHCCD161E %4.4X file[%d] no shadow file name\n"),
                dev->devnum, cckd->sfn+1);
        return -1;
    }

    /* Error if max number of shadow files exceeded */
    if (cckd->sfn+1 == CCKD_MAX_SF)
    {
        logmsg (_("HHCCD161E %4.4X file[%d] max shadow files exceeded\n"),
                dev->devnum, cckd->sfn+1);
        return -1;
    }

    /* Harden the current file */
    cckd_harden (dev);

    /* Open the new shadow file */
    if (cckd_open(dev, cckd->sfn+1, O_RDWR|O_CREAT|O_EXCL|O_BINARY,
                                      S_IRUSR | S_IWUSR | S_IRGRP) < 0)
        return -1;

    /* Read previous file's device header */
    if (cckd_read (dev, cckd->sfn, 0, &devhdr, CKDDASD_DEVHDR_SIZE) < 0)
        goto sf_new_error;

    /* Make sure identifier is CKD_S370 or FBA_S370 */
    devhdr.devid[4] = 'S';

    /* Write new file's device header */
    if (cckd_write (dev, cckd->sfn+1, 0, &devhdr, CKDDASD_DEVHDR_SIZE) < 0)
        goto sf_new_error;

    /* Build the compressed device header */
    memcpy (&cckd->cdevhdr[cckd->sfn+1], &cckd->cdevhdr[cckd->sfn], CCKDDASD_DEVHDR_SIZE);
    l1size = cckd->cdevhdr[cckd->sfn+1].numl1tab * CCKD_L1ENT_SIZE;
    cckd->cdevhdr[cckd->sfn+1].size =
    cckd->cdevhdr[cckd->sfn+1].used = CKDDASD_DEVHDR_SIZE + CCKDDASD_DEVHDR_SIZE + l1size;
    cckd->cdevhdr[cckd->sfn+1].free =
    cckd->cdevhdr[cckd->sfn+1].free_total =
    cckd->cdevhdr[cckd->sfn+1].free_largest =
    cckd->cdevhdr[cckd->sfn+1].free_number =
    cckd->cdevhdr[cckd->sfn+1].free_imbed = 0;

    /* Init the level 1 table */
    if ((cckd->l1[cckd->sfn+1] = cckd_malloc (dev, "l1", l1size)) == NULL)
        goto sf_new_error;
    memset (cckd->l1[cckd->sfn+1], 0xff, l1size);

    /* Make the new file active */
    cckd->sfn++;

    /* Harden the new file */
    if (cckd_harden (dev) < 0)
    {
        cckd->sfn--;
        goto sf_new_error;
    }

    /* Re-read the l1 to set l2bounds, l2ok */
    cckd_read_l1 (dev);

    return 0;

sf_new_error:
    cckd->l1[cckd->sfn+1] = cckd_free(dev, "l1", cckd->l1[cckd->sfn+1]);
    cckd_close (dev, cckd->sfn+1);
    cckd->open[cckd->sfn+1] = CCKD_OPEN_NONE;
    unlink (cckd_sf_name (dev, cckd->sfn+1));

    /* Re-read the l1 to set l2bounds, l2ok */
    cckd_read_l1 (dev);

    return -1;

} /* end function cckd_sf_new */

/*-------------------------------------------------------------------*/
/* Add a shadow file  (sf+)                                          */
/*-------------------------------------------------------------------*/
void *cckd_sf_add (void *data)
{
DEVBLK         *dev = data;             /* -> DEVBLK                 */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             syncio;                 /* Saved syncio bit          */

    if (dev == NULL)
    {
    int n = 0;
        for (dev=sysblk.firstdev; dev; dev=dev->nextdev)
            if (dev->cckd_ext)
            {
                logmsg( _("HHCCD207I Adding device %d:%4.4X\n"),
                          SSID_TO_LCSS(dev->ssid), dev->devnum );
                cckd_sf_add (dev);
                n++;
            }
        logmsg( _("HHCCD092I %d devices processed\n"), n );
        return NULL;
    }

    cckd = dev->cckd_ext;
    if (!cckd)
    {
        logmsg (_("HHCCD160E %4.4X not a cckd device\n"), dev->devnum);
        return NULL;
    }

    /* Disable synchronous I/O for the device */
    syncio = cckd_disable_syncio(dev);

    /* Schedule updated track entries to be written */
    obtain_lock (&cckd->iolock);
    if (cckd->merging)
    {
        dev->syncio = syncio;
        release_lock (&cckd->iolock);
        logmsg (_("HHCCD165W %4.4X error adding shadow file, "
                  "sf command busy on device\n"),
                dev->devnum,cckd->sfn);
        return NULL;
    }
    cckd->merging = 1;
    cckd_flush_cache (dev);
    while (cckd->wrpending || cckd->ioactive)
    {
        cckd->iowaiters++;
        wait_condition (&cckd->iocond, &cckd->iolock);
        cckd->iowaiters--;
        cckd_flush_cache (dev);
    }
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    release_lock (&cckd->iolock);

    /* Obtain control of the file */
    obtain_lock (&cckd->filelock);

    /* Harden the current file */
    cckd_harden (dev);

    /* Create a new shadow file */
    if (cckd_sf_new (dev) < 0) {
        logmsg (_("HHCCD161E %4.4X file[%d] error adding shadow file\n"),
                dev->devnum, cckd->sfn+1);
        goto cckd_sf_add_exit;
    }

    /* Re-open the previous file if opened read-write */
    if (cckd->open[cckd->sfn-1] == CCKD_OPEN_RW)
        cckd_open (dev, cckd->sfn-1, O_RDONLY|O_BINARY, 0);

    logmsg (_("HHCCD162I %4.4X file[%d] %s added\n"),
            dev->devnum, cckd->sfn, cckd_sf_name (dev, cckd->sfn));

cckd_sf_add_exit:

    /* Re-read the l1 to set l2bounds, l2ok */
    cckd_read_l1 (dev);

    release_lock (&cckd->filelock);

    obtain_lock (&cckd->iolock);
    cckd->merging = 0;
    if (cckd->iowaiters)
        broadcast_condition (&cckd->iocond);
    dev->syncio = syncio;
    release_lock (&cckd->iolock);

    cckd_sf_stats (dev);
    return NULL;
} /* end function cckd_sf_add */

/*-------------------------------------------------------------------*/
/* Remove a shadow file  (sf-)                                       */
/*-------------------------------------------------------------------*/
void *cckd_sf_remove (void *data)
{
DEVBLK         *dev = data;             /* -> DEVBLK                 */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             syncio;                 /* Saved syncio bit          */
int             rc;                     /* Return code               */
int             from_sfx, to_sfx;       /* From/to file index        */
int             fix;                    /* nullfmt index             */
int             add = 0;                /* 1=Add shadow file back    */
int             l2updated = 0;          /* 1=L2 table was updated    */
int             i,j;                    /* Loop indexes              */
int             merge, force;           /* Flags                     */
off_t           pos;                    /* File offset               */
size_t          len;                    /* Length to read/write      */
int             size;                   /* Image size                */
int             trk = -1;               /* Track being read/written  */
CCKD_L2ENT      from_l2[256],           /* Level 2 tables            */
                to_l2[256];
CCKD_L2ENT      new_l2;                 /* New level 2 table entry   */
BYTE            buf[65536];             /* Buffer                    */

    if (dev == NULL)
    {
    int n = 0;
        merge = cckdblk.sfmerge;
        force = cckdblk.sfforce;
        cckdblk.sfmerge = cckdblk.sfforce = 0;
        for (dev=sysblk.firstdev; dev; dev=dev->nextdev)
            if ((cckd = dev->cckd_ext))
            {
                logmsg( _("HHCCD179I Merging device %d:%4.4X\n"),
                          SSID_TO_LCSS(dev->ssid), dev->devnum );
                cckd->sfmerge = merge;
                cckd->sfforce = force;
                cckd_sf_remove (dev);
                n++;
            }
        logmsg( _("HHCCD092I %d devices processed\n"), n );
        return NULL;
    }

    cckd = dev->cckd_ext;
    if (!cckd)
    {
        logmsg (_("HHCCD170E %4.4X not a cckd device\n"),
                dev ? dev->devnum : 0);
        return NULL;
    }

    /* Set flags */
    merge = cckd->sfmerge || cckd->sfforce;
    force = cckd->sfforce;
    cckd->sfmerge = cckd->sfforce = 0;

    cckd_trace (dev, "merge starting: %s %s\n",
                merge ? "merge" : "nomerge", force ? "force" : "");

    /* Disable synchronous I/O for the device */
    syncio = cckd_disable_syncio(dev);

    /* Schedule updated track entries to be written */
    obtain_lock (&cckd->iolock);
    if (cckd->merging)
    {
        dev->syncio = syncio;
        release_lock (&cckd->iolock);
        logmsg (_("HHCCD175W %4.4X file[%d] merge failed, "
                  "sf command busy on device\n"),
                dev->devnum,cckd->sfn);
        return NULL;
    }
    cckd->merging = 1;
    cckd_flush_cache (dev);
    while (cckd->wrpending || cckd->ioactive)
    {
        cckd->iowaiters++;
        wait_condition (&cckd->iocond, &cckd->iolock);
        cckd->iowaiters--;
        cckd_flush_cache (dev);
    }
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    release_lock (&cckd->iolock);

    obtain_lock (&cckd->filelock);

    if (cckd->sfn == 0)
    {
        dev->syncio = syncio;
        release_lock (&cckd->filelock);
        logmsg (_("HHCCD171E %4.4X file[%d] cannot remove base file\n"),
                dev->devnum,cckd->sfn);
        cckd->merging = 0;
        return NULL;
    }

    from_sfx = cckd->sfn;
    to_sfx = cckd->sfn - 1;
    fix = cckd->cdevhdr[to_sfx].nullfmt;

    /* Harden the `from' file */
    if (cckd_harden (dev) < 0)
    {
        logmsg (_("HHCCD174E %4.4X file[%d] not merged, "
                "file[%d] not hardened\n"),
                dev->devnum, from_sfx, from_sfx);
        goto sf_remove_exit;
    }

    /* Attempt to re-open the `to' file read-write */
    cckd_close (dev, to_sfx);
    if (to_sfx > 0 || !dev->ckdrdonly || force)
        cckd_open (dev, to_sfx, O_RDWR|O_BINARY, 1);
    if (cckd->fd[to_sfx] < 0)
    {
        /* `from' file can't be opened read-write */
        cckd_open (dev, to_sfx, O_RDONLY|O_BINARY, 0);
        if (merge)
        {
            logmsg (_("HHCCD172E %4.4X file[%d] not merged, "
                    "file[%d] cannot be opened read-write%s\n"),
                    dev->devnum, from_sfx, to_sfx,
                    to_sfx == 0 && dev->ckdrdonly && !force
                    ? ", try `force'" : "");
            goto sf_remove_exit;
        }
        else
           add = 1;
    }
    else
    {
        /* `from' file opened read-write */
        cckd->sfn = to_sfx;
        if (cckd_chkdsk (dev, 0) < 0)
        {
            cckd->sfn = from_sfx;
            logmsg (_("HHCCD173E %4.4X file[%d] not merged, "
                    "file[%d] check failed\n"),
                    dev->devnum, to_sfx, to_sfx);
            goto sf_remove_exit;
        }
    }

    cckd->sfn = to_sfx;

    /* Perform backwards merge */
    if (merge)
    {
        cckd_trace (dev, "merging to file[%d]\n", to_sfx);

        /* Make the target file the active file */
        cckd->sfn = to_sfx;
        cckd->cdevhdr[to_sfx].options |= (CCKD_OPENED | CCKD_ORDWR);

        /* Loop for each level 1 table entry */
        for (i = 0; i < cckd->cdevhdr[from_sfx].numl1tab; i++)
        {
            l2updated = 0;
            /* Continue if from L2 doesn't exist */
            if (cckd->l1[from_sfx][i] == 0xffffffff
             || (cckd->l1[from_sfx][i] == 0 && cckd->l1[to_sfx][i] == 0))
                continue;

            trk = i*256;

            /* Read `from' l2 table */
            if (cckd->l1[from_sfx][i] == 0)
                memset (&from_l2, 0, CCKD_L2TAB_SIZE);
            else
            {
                pos = (off_t)cckd->l1[from_sfx][i];
                if (cckd_read(dev, from_sfx, pos, &from_l2, CCKD_L2TAB_SIZE) < 0)
                    goto sf_merge_error;
            }

            /* Read `to' l2 table */
            if (cckd->l1[to_sfx][i] == 0)
                memset (&to_l2, 0, CCKD_L2TAB_SIZE);
            else if (cckd->l1[to_sfx][i] == 0xffffffff)
                memset (&to_l2, 0xff, CCKD_L2TAB_SIZE);
            else
            {
                pos = (off_t)cckd->l1[to_sfx][i];
                if (cckd_read(dev, to_sfx, pos, &to_l2, CCKD_L2TAB_SIZE) < 0)
                    goto sf_merge_error;
            }

            /* Loop for each level 2 table entry */
            for (j = 0; j < 256; j++)
            {
                trk = i*256 + j;
                /* Continue if from L2 entry doesn't exist */
                if (from_l2[j].pos == 0xffffffff
                 || (from_l2[j].pos == 0 && to_l2[j].pos == 0))
                    continue;

                /* Read the `from' track/blkgrp image */
                len = (int)from_l2[j].len;
                if (len > CKDDASD_NULLTRK_FMTMAX)
                {
                    pos = (off_t)from_l2[j].pos;
                    if (cckd_read (dev, from_sfx, pos, buf, len) < 0)
                        goto sf_merge_error;

                    /* Get space for the `to' track/blkgrp image */
                    size = len;
                    if ((pos = cckd_get_space (dev, &size, CCKD_SIZE_EXACT)) < 0)
                        goto sf_merge_error;

                    new_l2.pos = (U32)pos;
                    new_l2.len = (U16)len;
                    new_l2.size = (U16)size;

                    /* Write the `to' track/blkgrp image */
                    if (cckd_write(dev, to_sfx, pos, buf, len) < 0)
                        goto sf_merge_error;
                }
                else
                {
                    new_l2.pos = 0;
                    new_l2.len = new_l2.size = (U16)len;
                }

                /* Release space occupied by old `to' entry */
                cckd_rel_space (dev, (off_t)to_l2[j].pos, (int)to_l2[j].len,
                                                          (int)to_l2[j].size);

                /* Update `to' l2 table entry */
                l2updated = 1;
                to_l2[j].pos = new_l2.pos;
                to_l2[j].len = new_l2.len;
                to_l2[j].size = new_l2.size;
            } /* for each level 2 table entry */

            /* Update the `to' level 2 table */
            if (l2updated)
            {
                l2updated = 0;
                pos = (off_t)cckd->l1[to_sfx][i];
                if (memcmp (&to_l2, &empty_l2[fix], CCKD_L2TAB_SIZE) == 0)
                {
                    cckd_rel_space (dev, pos, CCKD_L2TAB_SIZE, CCKD_L2TAB_SIZE);
                    pos = 0;
                }
                else
                {
                    size = CCKD_L2TAB_SIZE;
                    if (pos == 0 || pos == (off_t)0xffffffff)
                        if ((pos = cckd_get_space (dev, &size, CCKD_L2SPACE)) < 0)
                            goto sf_merge_error;
                    if (cckd_write(dev, to_sfx, pos, &to_l2, CCKD_L2TAB_SIZE) < 0)
                        goto sf_merge_error;
                } /* `to' level 2 table not null */

                /* Update the level 1 table index */
                cckd->l1[to_sfx][i] = (U32)pos;

                /* Flush free space */
                cckd_flush_space (dev);

            } /* Update level 2 table */

        } /* For each level 1 table entry */

        /* Validate the merge */
        cckd_harden (dev);
        cckd_chkdsk (dev, 0);
        cckd_read_init (dev);

    } /* if merge */

    /* Remove the old file */
    cckd_close (dev, from_sfx);
    cckd->l1[from_sfx] = cckd_free (dev, "l1", cckd->l1[from_sfx]);
    memset (&cckd->cdevhdr[from_sfx], 0, CCKDDASD_DEVHDR_SIZE);
    rc = unlink (cckd_sf_name (dev, from_sfx));

    /* Add the file back if necessary */
    if (add) rc = cckd_sf_new (dev) ;

    logmsg (_("HHCCD181I %4.4X shadow file [%d] successfully %s\n"),
            dev->devnum, from_sfx, merge ? "merged" : add ? "re-added" : "removed");

sf_remove_exit:

    /* Re-read the l1 to set l2bounds, l2ok */
    cckd_read_l1 (dev);

    release_lock (&cckd->filelock);

    obtain_lock (&cckd->iolock);
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    cckd->merging = 0;
    if (cckd->iowaiters)
        broadcast_condition (&cckd->iocond);
    dev->syncio = syncio;
    cckd_trace (dev, "merge complete%s\n","");
    release_lock (&cckd->iolock);

    cckd_sf_stats (dev);
    return NULL;

sf_merge_error:

    if (trk < 0)
        logmsg (_("HHCCD180E %4.4X file[%d] not merged, error during merge\n"),
                dev->devnum, from_sfx);
    else
        logmsg (_("HHCCD180E %4.4X file[%d] not merged, error processing trk %d\n"),
                dev->devnum, from_sfx, trk);

    if (l2updated && cckd->l1[to_sfx][i] && cckd->l1[to_sfx][i] != 0xffffffff)
    {
        l2updated = 0;
        pos = (off_t)cckd->l1[to_sfx][i];
        cckd_write(dev, to_sfx, pos, &to_l2, CCKD_L2TAB_SIZE);
    }
    cckd_harden(dev);
    cckd_chkdsk (dev, 2);
    cckd->sfn = from_sfx;
    cckd_harden(dev);
    cckd_chkdsk (dev, 2);
    goto sf_remove_exit;

} /* end function cckd_sf_remove */

/*-------------------------------------------------------------------*/
/* Check and compress a shadow file  (sfc)                           */
/*-------------------------------------------------------------------*/
void *cckd_sf_comp (void *data)
{
DEVBLK         *dev = data;             /* -> DEVBLK                 */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             syncio;                 /* Saved syncio bit          */
int             rc;                     /* Return code               */

    if (dev == NULL)
    {
    int n = 0;
        for (dev=sysblk.firstdev; dev; dev=dev->nextdev)
            if (dev->cckd_ext)
            {
                logmsg( _("HHCCD207I Compressing device %d:%4.4X\n"),
                          SSID_TO_LCSS(dev->ssid), dev->devnum );
                cckd_sf_comp (dev);
                n++;
            }
        logmsg( _("HHCCD092I %d devices processed\n"), n );
        return NULL;
    }

    cckd = dev->cckd_ext;
    if (!cckd)
    {
        logmsg (_("HHCCD205W %4.4X device is not a cckd device\n"), dev->devnum);
        return NULL;
    }

    /* Disable synchronous I/O for the device */
    syncio = cckd_disable_syncio(dev);

    /* schedule updated track entries to be written */
    obtain_lock (&cckd->iolock);
    if (cckd->merging)
    {
        dev->syncio = syncio;
        release_lock (&cckd->iolock);
        logmsg (_("HHCCD206W %4.4X file[%d] compress failed, "
                  "sf command busy on device\n"),
                dev->devnum,cckd->sfn);
        return NULL;
    }
    cckd->merging = 1;
    cckd_flush_cache (dev);
    while (cckd->wrpending || cckd->ioactive)
    {
        cckd->iowaiters++;
        wait_condition (&cckd->iocond, &cckd->iolock);
        cckd->iowaiters--;
        cckd_flush_cache (dev);
    }
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    release_lock (&cckd->iolock);

    /* obtain control of the file */
    obtain_lock (&cckd->filelock);

    /* harden the current file */
    cckd_harden (dev);

    /* Call the compress function */
    rc = cckd_comp (dev);

    /* Perform initial read */
    rc = cckd_read_init (dev);

    release_lock (&cckd->filelock);

    obtain_lock (&cckd->iolock);
    cckd->merging = 0;
    if (cckd->iowaiters)
        broadcast_condition (&cckd->iocond);
    dev->syncio = syncio;
    release_lock (&cckd->iolock);

    /* Display the shadow file statistics */
    cckd_sf_stats (dev);

    return NULL;
} /* end function cckd_sf_comp */

/*-------------------------------------------------------------------*/
/* Check a shadow file  (sfk)                                        */
/*-------------------------------------------------------------------*/
void *cckd_sf_chk (void *data)
{
DEVBLK         *dev = data;             /* -> DEVBLK                 */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             syncio;                 /* Saved syncio bit          */
int             rc;                     /* Return code               */
int             level = 2;              /* Check level               */

    if (dev == NULL)
    {
    int n = 0;
        level = cckdblk.sflevel;
        cckdblk.sflevel = 0;
        for (dev=sysblk.firstdev; dev; dev=dev->nextdev)
            if ((cckd = dev->cckd_ext))
            {
                logmsg( _("HHCCD207I Checking device %d:%4.4X level %d\n"),
                          SSID_TO_LCSS(dev->ssid), dev->devnum, level );
                cckd->sflevel = level;
                cckd_sf_chk (dev);
                n++;
            }
        logmsg( _("HHCCD092I %d devices processed\n"), n );
        return NULL;
    }

    cckd = dev->cckd_ext;
    if (!cckd)
    {
        logmsg (_("HHCCD205W %4.4X device is not a cckd device\n"), dev->devnum);
        return NULL;
    }

    level = cckd->sflevel;
    cckd->sflevel = 0;

    /* Disable synchronous I/O for the device */
    syncio = cckd_disable_syncio(dev);

    /* schedule updated track entries to be written */
    obtain_lock (&cckd->iolock);
    if (cckd->merging)
    {
        dev->syncio = syncio;
        release_lock (&cckd->iolock);
        logmsg (_("HHCCD206W %4.4X file[%d] check failed, "
                  "sf command busy on device\n"),
                dev->devnum,cckd->sfn);
        return NULL;
    }
    cckd->merging = 1;
    cckd_flush_cache (dev);
    while (cckd->wrpending || cckd->ioactive)
    {
        cckd->iowaiters++;
        wait_condition (&cckd->iocond, &cckd->iolock);
        cckd->iowaiters--;
        cckd_flush_cache (dev);
    }
    cckd_purge_cache (dev); cckd_purge_l2 (dev);
    dev->bufcur = dev->cache = -1;
    release_lock (&cckd->iolock);

    /* obtain control of the file */
    obtain_lock (&cckd->filelock);

    /* harden the current file */
    cckd_harden (dev);

    /* Call the chkdsk function */
    rc = cckd_chkdsk (dev, level);

    /* Perform initial read */
    rc = cckd_read_init (dev);

    release_lock (&cckd->filelock);

    obtain_lock (&cckd->iolock);
    cckd->merging = 0;
    if (cckd->iowaiters)
        broadcast_condition (&cckd->iocond);
    dev->syncio = syncio;
    release_lock (&cckd->iolock);

    /* Display the shadow file statistics */
    cckd_sf_stats (dev);

    return NULL;
} /* end function cckd_sf_chk */

/*-------------------------------------------------------------------*/
/* Display shadow file statistics   (sfd)                            */
/*-------------------------------------------------------------------*/
void *cckd_sf_stats (void *data)
{
DEVBLK         *dev = data;             /* -> DEVBLK                 */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
struct stat     st;                     /* File information          */
int             i;                      /* Index                     */
int             rc;                     /* Return code               */
char           *ost[] = {"  ", "ro", "rd", "rw"};
unsigned long long size=0,free=0;       /* Total size, free space    */
int             freenbr=0;              /* Total number free spaces  */

    if (dev == NULL)
    {
    int n = 0;
        for (dev=sysblk.firstdev; dev; dev=dev->nextdev)
            if (dev->cckd_ext)
            {
                logmsg( _("HHCCD208I Displaying device %d:%4.4X\n"),
                          SSID_TO_LCSS(dev->ssid), dev->devnum );
                cckd_sf_stats (dev);
                n++;
            }
        logmsg( _("HHCCD092I %d devices processed\n"), n );
        return NULL;
    }

    cckd = dev->cckd_ext;
    if (!cckd)
    {
        logmsg (_("HHCCD209W %4.4X device is not a cckd device\n"));
        return NULL;
    }

//  obtain_lock (&cckd->filelock);

    /* Calculate totals */
    rc = fstat (cckd->fd[0], &st);
    for (i = 0; i <= cckd->sfn; i++)
    {
        if (!i) size = st.st_size;
        else size += cckd->cdevhdr[i].size;
        free += cckd->cdevhdr[i].free_total;
        freenbr += cckd->cdevhdr[i].free_number;
    }

    /* header */
    logmsg (_("HHCCD210I           size free  nbr st   reads  writes l2reads    hits switches\n"));
    if (cckd->readaheads || cckd->misses)
    logmsg (_("HHCCD211I                                                  readaheads   misses\n"));
    logmsg (_("HHCCD212I --------------------------------------------------------------------\n"));

    /* total statistics */
    logmsg (_("HHCCD213I [*] %10" I64_FMT "d %3" I64_FMT "d%% %4d    %7d %7d %7d %7d  %7d\n"),
            size, (free * 100) / size, freenbr,
            cckd->totreads, cckd->totwrites, cckd->totl2reads,
            cckd->cachehits, cckd->switches);
    if (cckd->readaheads || cckd->misses)
    logmsg (_("HHCCD214I                                                     %7d  %7d\n"),
            cckd->readaheads, cckd->misses);

    /* base file statistics */
    logmsg (_("HHCCD215I %s\n"), dev->filename);
    logmsg (_("HHCCD216I [0] %10" I64_FMT "d %3" I64_FMT "d%% %4d %s %7d %7d %7d\n"),
            (long long)st.st_size,
            (long long)((long long)((long long)cckd->cdevhdr[0].free_total * 100) / st.st_size),
            cckd->cdevhdr[0].free_number, ost[cckd->open[0]],
            cckd->reads[0], cckd->writes[0], cckd->l2reads[0]);

    if (dev->dasdsfn != NULL && CCKD_MAX_SF > 0)
        logmsg (_("HHCCD217I %s\n"), cckd_sf_name(dev, -1));

    /* shadow file statistics */
    for (i = 1; i <= cckd->sfn; i++)
    {
        logmsg (_("HHCCD218I [%d] %10" I64_FMT "d %3" I64_FMT "d%% %4d %s %7d %7d %7d\n"),
                i, (long long)cckd->cdevhdr[i].size,
                (long long)((long long)((long long)cckd->cdevhdr[i].free_total * 100) / cckd->cdevhdr[i].size),
                cckd->cdevhdr[i].free_number, ost[cckd->open[i]],
                cckd->reads[i], cckd->writes[i], cckd->l2reads[i]);
    }
//  release_lock (&cckd->filelock);
    return NULL;
} /* end function cckd_sf_stats */

/*-------------------------------------------------------------------*/
/* Disable synchronous I/O for a device                              */
/*-------------------------------------------------------------------*/
int cckd_disable_syncio(DEVBLK *dev)
{
    if (!dev->syncio) return 0;
    obtain_lock(&dev->lock);
    while (dev->syncio_active)
    {
        release_lock(&dev->lock);
        usleep(500);
        obtain_lock(&dev->lock);
    }
    dev->syncio = 0;
    release_lock(&dev->lock);
    cckd_trace (dev, "syncio disabled%s\n","");
    return 1;
}

/*-------------------------------------------------------------------*/
/* Lock/unlock the device chain                                      */
/*-------------------------------------------------------------------*/
void cckd_lock_devchain(int flag)
{
//struct timespec tm;
//struct timeval  now;
//int             timeout;

    obtain_lock(&cckdblk.devlock);
    while ((flag && cckdblk.devusers != 0)
        || (!flag && cckdblk.devusers < 0))
    {
//gettimeofday(&now,NULL);
//tm.tv_sec = now.tv_sec + 2;
//tm.tv_nsec = now.tv_usec * 1000;
        cckdblk.devwaiters++;
        wait_condition(&cckdblk.devcond, &cckdblk.devlock);
//timeout = timed_wait_condition(&cckdblk.devcond, &cckdblk.devlock, &tm);
//if (timeout) cckd_print_itrace();
        cckdblk.devwaiters--;
    }
    if (flag) cckdblk.devusers--;
    else cckdblk.devusers++;
    release_lock(&cckdblk.devlock);
}
void cckd_unlock_devchain()
{
    obtain_lock(&cckdblk.devlock);
    if (cckdblk.devusers < 0) cckdblk.devusers++;
    else cckdblk.devusers--;
    if (cckdblk.devusers == 0 && cckdblk.devwaiters)
        signal_condition(&cckdblk.devcond);
    release_lock(&cckdblk.devlock);
}

/*-------------------------------------------------------------------*/
/* Garbage Collection thread                                         */
/*-------------------------------------------------------------------*/
void *cckd_gcol()
{
int             gcol;                   /* Identifier                */
int             rc;                     /* Return code               */
DEVBLK         *dev;                    /* -> device block           */
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
long long       size, fsiz;             /* File size, free size      */
struct timeval  tv_now;                 /* Time-of-day (as timeval)  */
time_t          tt_now;                 /* Time-of-day (as time_t)   */
struct timespec tm;                     /* Time-of-day to wait       */
int             gc;                     /* Garbage collection state  */
int             gctab[5]= {             /* default gcol parameters   */
                           4096,        /* critical  50%   - 100%    */
                           2048,        /* severe    25%   -  50%    */
                           1024,        /* moderate  12.5% -  25%    */
                            512,        /* light      6.3% -  12.5%  */
                            256};       /* none       0%   -   6.3%  */
//char *gcstates[] = {"critical","severe","moderate","light","none"};

    obtain_lock (&cckdblk.gclock);
    gcol = ++cckdblk.gcs;

    /* Return without messages if too many already started */
    if (gcol > cckdblk.gcmax)
    {
        --cckdblk.gcs;
        release_lock (&cckdblk.gclock);
        return NULL;
    }

    if (!cckdblk.batch)
    {
        logmsg (_("HHCCD003I Garbage collector thread started: tid="TIDPAT", pid=%d \n"),
              thread_id(), getpid());
    }

    while (gcol <= cckdblk.gcmax)
    {
        cckd_lock_devchain(0);
        /* Perform collection on each device */
        for (dev = cckdblk.dev1st; dev; dev = cckd->devnext)
        {
            cckd = dev->cckd_ext;
            obtain_lock (&cckd->iolock);

            /* Bypass if merging or stopping */
            if (cckd->merging || cckd->stopping)
            {
                release_lock (&cckd->iolock);
                continue;
            }

            /* Bypass if not opened read-write */
            if (cckd->open[cckd->sfn] != CCKD_OPEN_RW)
            {
                release_lock (&cckd->iolock);
                continue;
            }

            /* Free newbuf if it hasn't been used */
            if (!cckd->ioactive && !cckd->bufused && cckd->newbuf)
                cckd->newbuf = cckd_free (dev, "newbuf", cckd->newbuf);
            cckd->bufused = 0;

            /* If OPENED bit not on then flush if updated */
            if (!(cckd->cdevhdr[cckd->sfn].options & CCKD_OPENED))
            {
                if (cckd->updated) cckd_flush_cache (dev);
                release_lock (&cckd->iolock);
                continue;
            }

            /* Determine garbage state */
            size = (long long)cckd->cdevhdr[cckd->sfn].size;
            fsiz = (long long)cckd->cdevhdr[cckd->sfn].free_total;
            if      (fsiz >= (size = size/2)) gc = 0;
            else if (fsiz >= (size = size/2)) gc = 1;
            else if (fsiz >= (size = size/2)) gc = 2;
            else if (fsiz >= (size = size/2)) gc = 3;
            else gc = 4;

            /* Adjust the state based on the number of free spaces */
            if (cckd->cdevhdr[cckd->sfn].free_number >  800 && gc > 0) gc--;
            if (cckd->cdevhdr[cckd->sfn].free_number > 1800 && gc > 0) gc--;
            if (cckd->cdevhdr[cckd->sfn].free_number > 3000)           gc = 0;

            /* Set the size */
            if (cckdblk.gcparm > 0) size = gctab[gc] << cckdblk.gcparm;
            else if (cckdblk.gcparm < 0) size = gctab[gc] >> abs(cckdblk.gcparm);
            else size = gctab[gc];
            if (size > cckd->cdevhdr[cckd->sfn].used >> 10)
                size = cckd->cdevhdr[cckd->sfn].used >> 10;
            if (size < 64) size = 64;

            release_lock (&cckd->iolock);

            /* Call the garbage collector */
            cckd_gc_percolate (dev, (unsigned int)size);

            /* Schedule any updated tracks to be written */
            obtain_lock (&cckd->iolock);
            cckd_flush_cache (dev);
            while (cckdblk.fsync && cckd->wrpending)
            {
                cckd->iowaiters++;
                wait_condition (&cckd->iocond, &cckd->iolock);
                cckd->iowaiters--;
            }
            release_lock (&cckd->iolock);

            /* Sync the file */
            if (cckdblk.fsync && cckd->lastsync + 10 <= tv_now.tv_sec)
            {
                obtain_lock (&cckd->filelock);
                rc = fdatasync (cckd->fd[cckd->sfn]);
                cckd->lastsync = tv_now.tv_sec;
                release_lock (&cckd->filelock);
            }

            /* Flush the free space */
            if (cckd->cdevhdr[cckd->sfn].free_number)
            {
                obtain_lock (&cckd->filelock);
                cckd_flush_space (dev);
                release_lock (&cckd->filelock);
            }

        } /* for each cckd device */
        cckd_unlock_devchain();

        /* wait a bit */
        gettimeofday (&tv_now, NULL);
        tm.tv_sec = tv_now.tv_sec + cckdblk.gcwait;
        tm.tv_nsec = tv_now.tv_usec * 1000;
        tt_now = tv_now.tv_sec + ((tv_now.tv_usec + 500000)/1000000);
        cckd_trace (dev, "gcol wait %d seconds at %s",
                    cckdblk.gcwait, ctime (&tt_now));
        timed_wait_condition (&cckdblk.gccond, &cckdblk.gclock, &tm);
    }

    if (!cckdblk.batch)
    logmsg (_("HHCCD013I Garbage collector thread stopping: tid="TIDPAT", pid=%d\n"),
            thread_id(), getpid());

    cckdblk.gcs--;
    if (!cckdblk.gcs) signal_condition (&cckdblk.termcond);
    release_lock (&cckdblk.gclock);
    return NULL;
} /* end thread cckd_gcol */

/*-------------------------------------------------------------------*/
/* Garbage Collection -- Percolate algorithm                         */
/*-------------------------------------------------------------------*/
int cckd_gc_percolate(DEVBLK *dev, unsigned int size)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             rc;                     /* Return code               */
size_t          moved = 0;              /* Space moved               */
int             after = 0, a;           /* New space after old       */
int             sfx;                    /* File index                */
int             i, j, l;                /* Indexes                   */
int             flags;                  /* Write trkimg flags        */
off_t           fpos, upos;             /* File offsets              */
size_t          flen, ulen, len;        /* Lengths                   */
int             trk;                    /* Track number              */
int             l1x,l2x;                /* Table Indexes             */
CCKD_L2ENT      l2;                     /* Copied level 2 entry      */
BYTE            buf[256*1024];          /* Buffer                    */

    cckd = dev->cckd_ext;
    size = size << 10;

    /* Debug */
    if (cckdblk.itracen)
    {
        cckd_trace (dev, "gcperc size %d 1st 0x%x nbr %d largest %u\n",
                    size, cckd->cdevhdr[cckd->sfn].free,
                    cckd->cdevhdr[cckd->sfn].free_number,
                    cckd->cdevhdr[cckd->sfn].free_largest);
        fpos = (off_t)cckd->cdevhdr[cckd->sfn].free;
        for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
        {
            cckd_trace (dev, "gcperc free[%4d]:%8.8x end %8.8x len %10d%cpend %d\n",
                        i,(int)fpos,(int)(fpos+cckd->free[i].len),(int)cckd->free[i].len,
                        fpos+(int)cckd->free[i].len == (int)cckd->free[i].pos ?
                                '*' : ' ',cckd->free[i].pending);
            fpos = cckd->free[i].pos;
        }
    }

    if (!cckd->l2ok)
        cckd_gc_l2(dev, buf);

    /* garbage collection cycle */
    while (moved < size && after < 4)
    {
        obtain_lock (&cckd->filelock);
        sfx = cckd->sfn;

        /* Exit if no more free space */
        if (cckd->cdevhdr[sfx].free_total == 0)
        {
            release_lock (&cckd->filelock);
            return moved;
        }

        /* Make sure the free space chain is built */
        if (!cckd->free) cckd_read_fsp (dev);

        /* Find a space to start with */
        l = -1;
        upos = ulen = flen = 0;
        fpos = cckd->cdevhdr[sfx].free;

        /* First non-pending free space */
        for (i = cckd->free1st; i >= 0; i = cckd->free[i].next)
        {
            if (!cckd->free[i].pending)
            {
                flen += cckd->free[i].len;
                break;
            }
            fpos = cckd->free[i].pos;
        }

        /* Continue to largest if non-zero `after' */
        for ( ; i >= 0 && after; i = cckd->free[i].next)
        {
            l = i;
            if (!cckd->free[i].pending) flen += cckd->free[i].len;
            if (cckd->free[i].len == cckd->cdevhdr[sfx].free_largest)
                break;
            fpos = cckd->free[i].pos;
        }

        /* Skip following free spaces */
        for ( ; i >= 0; i = cckd->free[i].next)
        {
            if (!cckd->free[i].pending) flen += cckd->free[i].len;
            if (fpos + cckd->free[i].len != cckd->free[i].pos) break;
            fpos = cckd->free[i].pos;
        }

        /* Space preceding largest if largest is at the end */
        if (i < 0 && l >= 0)
        {
            if (!cckd->free[l].pending) flen -= cckd->free[i].len;
            for (i = cckd->free[l].prev; i >= 0; i = cckd->free[i].prev)
            {
                fpos = cckd->free[i].prev >= 0
                     ? cckd->free[cckd->free[i].prev].pos
                     : cckd->cdevhdr[sfx].free;
                if (fpos + cckd->free[i].len < cckd->free[i].pos) break;
                if (!cckd->free[i].pending) flen -= cckd->free[i].len;
            }
        }

        /* Calculate the offset/length of the used space.
         * If only imbedded free space is left, then start
         * with the first used space that is not an l2 table.
         */
        if (i >= 0)
        {
            upos = fpos + cckd->free[i].len;
            ulen = (cckd->free[i].pos ? cckd->free[i].pos : cckd->cdevhdr[sfx].size) - upos;
        }
        else if (!cckd->cdevhdr[sfx].free_number && cckd->cdevhdr[sfx].free_imbed)
        {
            upos = (off_t)(CCKD_L1TAB_POS + cckd->cdevhdr[sfx].numl1tab * CCKD_L1ENT_SIZE);
            while (1)
            {
                for (i = 0; i < cckd->cdevhdr[sfx].numl1tab; i++)
                    if (cckd->l1[sfx][i] == (U32)upos)
                       break;
                if (i >= cckd->cdevhdr[sfx].numl1tab)
                    break;
                upos += CCKD_L2TAB_SIZE;
            }
            ulen = cckd->cdevhdr[sfx].size - upos;
        }

        /* Return if no applicable used space */
        if (ulen == 0)
        {
            cckd_trace (dev, "gcperc no applicable space, moved %ld\n",(long)moved);
            release_lock (&cckd->filelock);
            return moved;
        }

        /* Reduce ulen size to minimize `after' relocations */
        if (ulen > flen + 65536) ulen = flen + 65536;
        if (ulen > sizeof(buf))  ulen = sizeof(buf);

        cckd_trace (dev, "gcperc selected space 0x%" I64_FMT "x len %ld\n",
                    (long long)upos, (long)ulen);

        if (cckd_read (dev, sfx, upos, buf, ulen) < 0)
            goto cckd_gc_perc_error;

        /* Process each space in the buffer */
        flags = cckd->cdevhdr[sfx].free_number < 100 ? CCKD_SIZE_EXACT : CCKD_SIZE_ANY;
        for (i = a = 0; i + CKDDASD_TRKHDR_SIZE <= (int)ulen; i += len)
        {
            /* Check for level 2 table */
            for (j = 0; j < cckd->cdevhdr[sfx].numl1tab; j++)
                if (cckd->l1[sfx][j] == (U32)(upos + i)) break;

            if (j < cckd->cdevhdr[sfx].numl1tab)
            {
                /* Moving a level 2 table */
                len = CCKD_L2TAB_SIZE;
                if (i + len > ulen) break;
                cckd_trace (dev, "gcperc move l2tab[%d] at pos 0x%" I64_FMT "x len %ld\n",
                            j, (unsigned long long)(upos + i), (long)len);

                /* Make the level 2 table active */
                if (cckd_read_l2 (dev, sfx, j) < 0)
                    goto cckd_gc_perc_error;

                /* Write the level 2 table */
                if (cckd_write_l2 (dev) < 0)
                    goto cckd_gc_perc_error;
            }
            else
            {
                /* Moving a track image */
                if ((trk = cckd_cchh (dev, buf + i, -1)) < 0)
                    goto cckd_gc_perc_space_error;

                l1x = trk >> 8;
                l2x = trk & 0xff;

                /* Read the lookup entry for the track */
                if (cckd_read_l2ent (dev, &l2, trk) < 0)
                    goto cckd_gc_perc_error;
                if (l2.pos != (U32)(upos + i))
                    goto cckd_gc_perc_space_error;
                len = (int)l2.size;
                if (i + l2.len > (int)ulen) break;

                cckd_trace (dev, "gcperc move trk %d at pos 0x%" I64_FMT "x len %d\n",
                            trk, (long long)(upos + i), (int)l2.len);

                /* Relocate the track image somewhere else */
                if ((rc = cckd_write_trkimg (dev, buf + i, (int)l2.len, trk, flags)) < 0)
                    goto cckd_gc_perc_error;
                a += rc;
            }
        } /* for each space in the used space */

        /* Set `after' to 1 if first time space was relocated after */
        after += after ? a : (a > 0);
        moved += i;

        cckdblk.stats_gcolmoves++;
        cckdblk.stats_gcolbytes += i;

        release_lock (&cckd->filelock);

    } /* while (moved < size) */

    cckd_trace (dev, "gcperc moved %ld 1st 0x%x nbr %d\n", (long)moved,
                cckd->cdevhdr[cckd->sfn].free,cckd->cdevhdr[cckd->sfn].free_number);
    return moved;

cckd_gc_perc_space_error:

    logmsg (_("HHCCD190E %4.4X file[%d] offset 0x%" I64_FMT "x unknown space: "
              "%2.2x%2.2x%2.2x%2.2x%2.2x\n"),
            dev->devnum,cckd->sfn,(long long)(upos + i),
            buf[i], buf[i+1],buf[i+2], buf[i+3], buf[i+4]);
    cckd->cdevhdr[cckd->sfn].options |= CCKD_SPERRS;
    cckd_print_itrace();

cckd_gc_perc_error:

    cckd_trace (dev, "gcperc exiting due to error, moved %ld\n", (long)moved);
    release_lock (&cckd->filelock);
    return moved;

} /* end function cckd_gc_percolate */

/*-------------------------------------------------------------------*/
/* Garbage Collection -- Reposition level 2 tables                   */
/*                                                                   */
/* General idea is to relocate all level 2 tables as close to the    */
/* beginning of the file as possible.  This can help speed up, for   */
/* example, chkdsk processing.                                       */
/*                                                                   */
/* If any level 2 tables reside outside of the bounds (that is, if   */
/* any level 2 table could be moved closer to the beginning of the   */
/* file) then first we relocate all track images within the bounds.  */
/* Note that cckd_get_space will not allocate space within the       */
/* the bounds for track images.  Next we try to relocate all level 2 */
/* tables outside the bounds.  This may take a few iterations for    */
/* the freed space within the bounds to become non-pending.          */
/*                                                                   */
/* The bounds can change as level 2 tables are added or removed.     */
/* cckd_read_l1 sets the bounds and they are adjusted by             */
/* cckd_write_l2.                                                    */
/*-------------------------------------------------------------------*/
int cckd_gc_l2(DEVBLK *dev, BYTE *buf)
{
CCKDDASD_EXT   *cckd;                   /* -> cckd extension         */
int             sfx;                    /* Shadow file index         */
int             i, j;                   /* Work variables            */
int             trk;                    /* Track number              */
int             len;                    /* Track length              */
off_t           pos, fpos;              /* File offsets              */

    cckd = dev->cckd_ext;

    obtain_lock (&cckd->filelock);
    sfx = cckd->sfn;

    if (cckd->l2ok || cckd->cdevhdr[cckd->sfn].free_total == 0)
        goto cckd_gc_l2_exit;

    /* Find any level 2 table out of bounds */
    for (i = 0; i < cckd->cdevhdr[sfx].numl1tab; i++)
        if (cckd->l1[sfx][i] != 0 && cckd->l1[sfx][i] != 0xffffffff
         && cckd->l2bounds - CCKD_L2TAB_SIZE < (off_t)cckd->l1[sfx][i])
            break;

    /* Return OK if no l2 tables out of bounds */
    if (i >= cckd->cdevhdr[sfx].numl1tab)
        goto cckd_gc_l2_exit_ok;

    /* Relocate all track images within the bounds */
    pos = CCKD_L1TAB_POS + (cckd->cdevhdr[sfx].numl1tab * CCKD_L1ENT_SIZE);
    i = cckd->free1st;
    fpos = (off_t)cckd->cdevhdr[sfx].free;
    while (pos < cckd->l2bounds)
    {
        if (i >= 0 && pos == fpos)
        {
            pos += cckd->free[i].len;
            fpos = (off_t)cckd->free[i].pos;
            i = cckd->free[i].next;
            j = 0;
        }
        else
        {
            for (j = 0; j < cckd->cdevhdr[sfx].numl1tab; j++)
                if (pos == (off_t)cckd->l1[sfx][j])
                {
                    pos += CCKD_L2TAB_SIZE;
                    break;
                }
        }
        if (j >= cckd->cdevhdr[sfx].numl1tab)
        {
            /* Found a track to relocate */
            if (cckd_read (dev, sfx, pos, buf, CKDDASD_TRKHDR_SIZE) < 0)
                goto cckd_gc_l2_exit;
            if ((trk = cckd_cchh (dev, buf, -1)) < 0)
                goto cckd_gc_l2_exit;
            cckd_trace (dev, "gc_l2 relocate trk[%d] offset 0x%x\n", trk, pos);
            if ((len = cckd_read_trkimg (dev, buf, trk, NULL)) < 0)
               goto cckd_gc_l2_exit;
            if (cckd_write_trkimg (dev, buf, len, trk, CCKD_SIZE_EXACT) < 0)
               goto cckd_gc_l2_exit;
            /* Start over */
            pos = CCKD_L1TAB_POS + (cckd->cdevhdr[sfx].numl1tab * CCKD_L1ENT_SIZE);
            i = cckd->free1st;
            fpos = (off_t)cckd->cdevhdr[sfx].free;
        }
    }

    do {
        /* Find a level 2 table to relocate */
        i = cckd->free1st;
        fpos = (off_t)cckd->cdevhdr[sfx].free;
        cckd_trace (dev, "gc_l2 first free[%d] pos 0x%x len %d pending %d\n",
                    i, (int)fpos, i >= 0 ? (int)cckd->free[i].len : -1,
                    i >= 0 ? cckd->free[i].pending : -1);
        if (i < 0 || fpos >= cckd->l2bounds || cckd->free[i].pending)
            goto cckd_gc_l2_exit;

        if (cckd->free[i].len < CCKD_L2TAB_SIZE
         || (cckd->free[i].len != CCKD_L2TAB_SIZE
          && cckd->free[i].len < CCKD_L2TAB_SIZE + CCKD_FREEBLK_SIZE
            )
           )
        {
            for (i = 0; i < cckd->cdevhdr[sfx].numl1tab; i++)
                if (fpos + cckd->free[i].len == (off_t)cckd->l1[sfx][i])
                    break;
        }
        else
        {
            for (i = 0; i < cckd->cdevhdr[sfx].numl1tab; i++)
                if (cckd->l2bounds - CCKD_L2TAB_SIZE < (off_t)cckd->l1[sfx][i]
                 && cckd->l1[sfx][i] != 0xffffffff)
                    break;
        }

        if (i < cckd->cdevhdr[sfx].numl1tab)
        {
            cckd_trace (dev, "gc_l2 relocate l2[%d] pos 0x%x\n",
                        i, cckd->l1[sfx][i]);
            if (cckd_read_l2 (dev, sfx, i) < 0)
                goto cckd_gc_l2_exit;
            if (cckd_write_l2 (dev) < 0)
                goto cckd_gc_l2_exit;
        }
    } while (i < cckd->cdevhdr[sfx].numl1tab);

cckd_gc_l2_exit:
    release_lock (&cckd->filelock);
    return 0;

cckd_gc_l2_exit_ok:
    cckd_trace (dev, "gc_l2 ok%s\n", "");
    cckd->l2ok = 1;
    goto cckd_gc_l2_exit;
}

/*-------------------------------------------------------------------*/
/* Find device by devnum                                             */
/*-------------------------------------------------------------------*/
DEVBLK *cckd_find_device_by_devnum (U16 devnum)
{
DEVBLK       *dev;
CCKDDASD_EXT *cckd;

    cckd_lock_devchain (0);
    for (dev = cckdblk.dev1st; dev; dev = cckd->devnext)
    {
        if (dev->devnum == devnum) break;
        cckd = dev->cckd_ext;
    }
    cckd_unlock_devchain ();
    return dev;
} /* end function cckd_find_device_by_devnum */

/*-------------------------------------------------------------------*/
/* Uncompress a track image                                          */
/*-------------------------------------------------------------------*/
BYTE *cckd_uncompress (DEVBLK *dev, BYTE *from, int len, int maxlen,
                       int trk)
{
CCKDDASD_EXT   *cckd;
BYTE           *to = NULL;                /* Uncompressed buffer     */
int             newlen;                   /* Uncompressed length     */
BYTE            comp;                     /* Compression type        */
static char    *compress[] = {"none", "zlib", "bzip2"};

    cckd = dev->cckd_ext;

    cckd_trace (dev, "uncompress comp %d len %d maxlen %d trk %d\n",
                from[0] & CCKD_COMPRESS_MASK, len, maxlen, trk);

    /* Extract compression type */
    comp = (from[0] & CCKD_COMPRESS_MASK);

    /* Get a buffer to uncompress into */
    if (comp != CCKD_COMPRESS_NONE && cckd->newbuf == NULL)
    {
        cckd->newbuf = cckd_malloc (dev, "newbuf", maxlen);
        if (cckd->newbuf == NULL)
            return NULL;
    }

    /* Uncompress the track image */
    switch (comp) {

    case CCKD_COMPRESS_NONE:
        newlen = cckd_trklen (dev, from);
        to = from;
        break;
    case CCKD_COMPRESS_ZLIB:
        to = cckd->newbuf;
        newlen = cckd_uncompress_zlib (dev, to, from, len, maxlen);
        break;
    case CCKD_COMPRESS_BZIP2:
        to = cckd->newbuf;
        newlen = cckd_uncompress_bzip2 (dev, to, from, len, maxlen);
        break;
    default:
        newlen = -1;
        break;
    }

    /* Validate the uncompressed track image */
    newlen = cckd_validate (dev, to, trk, newlen);

    /* Return if successful */
    if (newlen > 0)
    {
        if (to != from)
        {
            cckd->newbuf = from;
            cckd->bufused = 1;
        }
        return to;
    }

    /* Get a buffer now if we haven't gotten one */
    if (cckd->newbuf == NULL)
    {
        cckd->newbuf = cckd_malloc (dev, "newbuf2", maxlen);
        if (cckd->newbuf == NULL)
            return NULL;
    }

    /* Try each uncompression routine in turn */

    /* uncompressed */
    newlen = cckd_trklen (dev, from);
    newlen = cckd_validate (dev, from, trk, newlen);
    if (newlen > 0)
        return from;

    /* zlib compression */
    to = cckd->newbuf;
    newlen = cckd_uncompress_zlib (dev, to, from, len, maxlen);
    newlen = cckd_validate (dev, to, trk, newlen);
    if (newlen > 0)
    {
        cckd->newbuf = from;
        cckd->bufused = 1;
        return to;
    }

    /* bzip2 compression */
    to = cckd->newbuf;
    newlen = cckd_uncompress_bzip2 (dev, to, from, len, maxlen);
    newlen = cckd_validate (dev, to, trk, newlen);
    if (newlen > 0)
    {
        cckd->newbuf = from;
        cckd->bufused = 1;
        return to;
    }

    /* Unable to uncompress */
    logmsg (_("HHCCD193E %4.4X file[%d] uncompress error trk %d: %2.2x%2.2x%2.2x%2.2x%2.2x\n"),
            dev->devnum, cckd->sfn, trk, from[0], from[1], from[2], from[3], from[4]);
    if (comp & ~cckdblk.comps)
        logmsg (_("HHCCD194E %4.4X file[%d] %s compression not supported\n"),
                dev->devnum, cckd->sfn, compress[comp]);
    return NULL;
}

int cckd_uncompress_zlib (DEVBLK *dev, BYTE *to, BYTE *from, int len, int maxlen)
{
#if defined(HAVE_LIBZ)
unsigned long newlen;
int rc;

    UNREFERENCED(dev);
    memcpy (to, from, CKDDASD_TRKHDR_SIZE);
    newlen = maxlen - CKDDASD_TRKHDR_SIZE;
    rc = uncompress(&to[CKDDASD_TRKHDR_SIZE], &newlen,
                &from[CKDDASD_TRKHDR_SIZE], len - CKDDASD_TRKHDR_SIZE);
    if (rc == Z_OK)
    {
        newlen += CKDDASD_TRKHDR_SIZE;
        to[0] = 0;
    }
    else
        newlen = -1;

    cckd_trace (dev, "uncompress zlib newlen %d rc %d\n",(int)newlen,rc);

    return (int)newlen;
#else
    UNREFERENCED(dev);
    UNREFERENCED(to);
    UNREFERENCED(from);
    UNREFERENCED(len);
    UNREFERENCED(maxlen);
    return -1;
#endif
}
int cckd_uncompress_bzip2 (DEVBLK *dev, BYTE *to, BYTE *from, int len, int maxlen)
{
#if defined(CCKD_BZIP2)
unsigned int newlen;
int rc;

    UNREFERENCED(dev);
    memcpy (to, from, CKDDASD_TRKHDR_SIZE);
    newlen = maxlen - CKDDASD_TRKHDR_SIZE;
    rc = BZ2_bzBuffToBuffDecompress (
                (void *)&to[CKDDASD_TRKHDR_SIZE], &newlen,
                (void *)&from[CKDDASD_TRKHDR_SIZE], len - CKDDASD_TRKHDR_SIZE,
                0, 0);
    if (rc == BZ_OK)
    {
        newlen += CKDDASD_TRKHDR_SIZE;
        to[0] = 0;
    }
    else
        newlen = -1;

    cckd_trace (dev, "uncompress bz2 newlen %d rc %d\n",newlen,rc);

    return (int)newlen;
#else
    UNREFERENCED(dev);
    UNREFERENCED(to);
    UNREFERENCED(from);
    UNREFERENCED(len);
    UNREFERENCED(maxlen);
    return -1;
#endif
}

/*-------------------------------------------------------------------*/
/* Compress a track image                                            */
/*-------------------------------------------------------------------*/
int cckd_compress (DEVBLK *dev, BYTE **to, BYTE *from, int len,
                   int comp, int parm)
{
int newlen;

    switch (comp) {
    case CCKD_COMPRESS_NONE:
        newlen = cckd_compress_none (dev, to, from, len, parm);
        break;
    case CCKD_COMPRESS_ZLIB:
        newlen = cckd_compress_zlib (dev, to, from, len, parm);
        break;
    case CCKD_COMPRESS_BZIP2:
        newlen = cckd_compress_bzip2 (dev, to, from, len, parm);
        break;
    default:
        newlen = cckd_compress_bzip2 (dev, to, from, len, parm);
        break;
    }
    return newlen;
}
int cckd_compress_none (DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm)
{
    UNREFERENCED(dev);
    UNREFERENCED(parm);
    *to = from;
    from[0] = CCKD_COMPRESS_NONE;
    return len;
}
int cckd_compress_zlib (DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm)
{
#if defined(HAVE_LIBZ)
unsigned long newlen;
int rc;
BYTE *buf;

    UNREFERENCED(dev);
    buf = *to;
    from[0] = CCKD_COMPRESS_NONE;
    memcpy (buf, from, CKDDASD_TRKHDR_SIZE);
    buf[0] = CCKD_COMPRESS_ZLIB;
    newlen = 65535 - CKDDASD_TRKHDR_SIZE;
    rc = compress2 (&buf[CKDDASD_TRKHDR_SIZE], &newlen,
                    &from[CKDDASD_TRKHDR_SIZE], len - CKDDASD_TRKHDR_SIZE,
                    parm);
    newlen += CKDDASD_TRKHDR_SIZE;
    if (rc != Z_OK || (int)newlen >= len)
    {
        *to = from;
        newlen = len;
    }
    return (int)newlen;
#else

#if defined(CCKD_BZIP2)
    return cckd_compress_bzip2 (dev, to, from, len, parm);
#else
    return cckd_compress_none (dev, to, from, len, parm);
#endif

#endif
}
int cckd_compress_bzip2 (DEVBLK *dev, BYTE **to, BYTE *from, int len, int parm)
{
#if defined(CCKD_BZIP2)
unsigned int newlen;
int rc;
BYTE *buf;

    UNREFERENCED(dev);
    buf = *to;
    from[0] = CCKD_COMPRESS_NONE;
    memcpy (buf, from, CKDDASD_TRKHDR_SIZE);
    buf[0] = CCKD_COMPRESS_BZIP2;
    newlen = 65535 - CKDDASD_TRKHDR_SIZE;
    rc = BZ2_bzBuffToBuffCompress (
                    (void *)&buf[CKDDASD_TRKHDR_SIZE], &newlen,
                    (void *)&from[CKDDASD_TRKHDR_SIZE], len - CKDDASD_TRKHDR_SIZE,
                    parm >= 1 && parm <= 9 ? parm : 5, 0, 0);
    newlen += CKDDASD_TRKHDR_SIZE;
    if (rc != BZ_OK || newlen >= (unsigned int)len)
    {
        *to = from;
        newlen = len;
    }
    return newlen;
#else
    return cckd_compress_zlib (dev, to, from, len, parm);
#endif
}

/*-------------------------------------------------------------------*/
/* cckd command help                                                 */
/*-------------------------------------------------------------------*/
void cckd_command_help()
{
    logmsg ("cckd command parameters:\n"
             "help\t\tDisplay help message\n"
             "stats\t\tDisplay cckd statistics\n"
             "opts\t\tDisplay cckd options\n"
             "comp=<n>\t\tOverride compression\t\t(-1 .. 2)\n"
             "compparm=<n>\tOverride compression parm\t\t(-1 .. 9)\n"
             "ra=<n>\t\tSet number readahead threads\t\t(1 .. 9)\n"
             "raq=<n>\t\tSet readahead queue size\t\t(0 .. 16)\n"
             "rat=<n>\t\tSet number tracks to read ahead\t\t(0 .. 16)\n"
             "wr=<n>\t\tSet number writer threads\t\t(1 .. 9)\n"
             "gcint=<n>\tSet garbage collector interval (sec)\t(1 .. 60)\n"
             "gcparm=<n>\tSet garbage collector parameter\t\t(-8 .. 8)\n"
             "\t\t    (least agressive ... most aggressive)\n"
             "nostress=<n>\t1=Disable stress writes\n"
             "freepend=<n>\tSet free pending cycles\t\t\t(-1 .. 4)\n"
             "fsync=<n>\t1=Enable fsync()\n"
             "trace=<n>\tSet trace table size\t\t\t(0 .. 200000)\n"
            );
} /* end function cckd_command_help */

/*-------------------------------------------------------------------*/
/* cckd command stats                                                */
/*-------------------------------------------------------------------*/
void cckd_command_opts()
{
    logmsg ("comp=%d,compparm=%d,ra=%d,raq=%d,rat=%d,"
             "wr=%d,gcint=%d,gcparm=%d,nostress=%d,\n"
             "\tfreepend=%d,fsync=%d,trace=%d,linuxnull=%d\n",
             cckdblk.comp == 0xff ? -1 : cckdblk.comp,
             cckdblk.compparm, cckdblk.ramax,
             cckdblk.ranbr, cckdblk.readaheads,
             cckdblk.wrmax, cckdblk.gcwait,
             cckdblk.gcparm, cckdblk.nostress, cckdblk.freepend,
             cckdblk.fsync, cckdblk.itracen, cckdblk.linuxnull);
} /* end function cckd_command_opts */

/*-------------------------------------------------------------------*/
/* cckd command stats                                                */
/*-------------------------------------------------------------------*/
void cckd_command_stats()
{
    logmsg("reads....%10" I64_FMT "d Kbytes...%10" I64_FMT "d writes...%10" I64_FMT "d Kbytes...%10" I64_FMT "d\n"
            "readaheads%9" I64_FMT "d misses...%10" I64_FMT "d syncios..%10" I64_FMT "d misses...%10" I64_FMT "d\n"
            "switches.%10" I64_FMT "d l2 reads.%10" I64_FMT "d              stress writes...%10" I64_FMT "d\n"
            "cachehits%10" I64_FMT "d misses...%10" I64_FMT "d l2 hits..%10" I64_FMT "d misses...%10" I64_FMT "d\n"
            "waits                                   i/o......%10" I64_FMT "d cache....%10" I64_FMT "d\n"
            "garbage collector   moves....%10" I64_FMT "d Kbytes...%10" I64_FMT "d\n",
            cckdblk.stats_reads, cckdblk.stats_readbytes >> 10,
            cckdblk.stats_writes, cckdblk.stats_writebytes >> 10,
            cckdblk.stats_readaheads, cckdblk.stats_readaheadmisses,
            cckdblk.stats_syncios, cckdblk.stats_synciomisses,
            cckdblk.stats_switches, cckdblk.stats_l2reads,
            cckdblk.stats_stresswrites,
            cckdblk.stats_cachehits, cckdblk.stats_cachemisses,
            cckdblk.stats_l2cachehits, cckdblk.stats_l2cachemisses,
            cckdblk.stats_iowaits, cckdblk.stats_cachewaits,
            cckdblk.stats_gcolmoves, cckdblk.stats_gcolbytes >> 10);
} /* end function cckd_command_stats */

/*-------------------------------------------------------------------*/
/* cckd command debug                                                */
/*-------------------------------------------------------------------*/
void cckd_command_debug()
{
}

/*-------------------------------------------------------------------*/
/* cckd command processor                                            */
/*-------------------------------------------------------------------*/
DLL_EXPORT int cckd_command(char *op, int cmd)
{
char  *kw, *p, c = '\0', buf[256];
int   val, opts = 0;

    /* Display help for null operand */
    if (op == NULL)
    {
        if (memcmp (&cckdblk.id, "CCKDBLK ", sizeof(cckdblk.id)) == 0 && cmd)
            cckd_command_help();
        return 0;
    }

    strcpy(buf, op);
    op = buf;

    /* Initialize the global cckd block if necessary */
    if (memcmp (&cckdblk.id, "CCKDBLK ", sizeof(cckdblk.id)))
        cckddasd_init (0, NULL);

    while (op)
    {
        /* Operands are delimited by commas */
        kw = op;
        op = strchr (op, ',');
        if (op) *op++ = '\0';

        /* Check for keyword = value */
        if ((p = strchr (kw, '=')))
        {
            *p++ = '\0';
            sscanf (p, "%d%c", &val, &c);
        }
        else val = -77;

        /* Parse the keyword */
        if (strcasecmp (kw, "help") == 0)
        {
            if (!cmd) return 0;
            cckd_command_help();
        }
        else if (strcasecmp (kw, "stats") == 0)
        {
            if (!cmd) return 0;
            cckd_command_stats ();
        }
        else if (strcasecmp (kw, "opts") == 0)
        {
            if (!cmd) return 0;
            cckd_command_opts();
        }
        else if (strcasecmp (kw, "debug") == 0)
        {
            if (!cmd) return 0;
            cckd_command_debug();
        }
        else if (strcasecmp (kw, "comp") == 0)
        {
            if (val < -1 || (val & ~cckdblk.comps) || c != '\0')
            {
                logmsg ("Invalid value for comp=\n");
                return -1;
            }
            else
            {
                cckdblk.comp = val < 0 ? 0xff : val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "compparm") == 0)
        {
            if (val < -1 || val > 9 || c != '\0')
            {
                logmsg ("Invalid value for compparm=\n");
                return -1;
            }
            else
            {
                cckdblk.compparm = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "ra") == 0)
        {
            if (val < CCKD_MIN_RA || val > CCKD_MAX_RA || c != '\0')
            {
                logmsg ("Invalid value for ra=\n");
                return -1;
            }
            else
            {
                cckdblk.ramax = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "raq") == 0)
        {
            if (val < 0 || val > CCKD_MAX_RA_SIZE || c != '\0')
            {
                logmsg ("Invalid value for raq=\n");
                return -1;
            }
            else
            {
                cckdblk.ranbr = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "rat") == 0)
        {
            if (val < 0 || val > CCKD_MAX_RA_SIZE || c != '\0')
            {
                logmsg ("Invalid value for rat=\n");
                return -1;
            }
            else
            {
                cckdblk.readaheads = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "wr") == 0)
        {
            if (val < CCKD_MIN_WRITER || val > CCKD_MAX_WRITER || c != '\0')
            {
                logmsg ("Invalid value for wr=\n");
                return -1;
            }
            else
            {
                cckdblk.wrmax = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "gcint") == 0)
        {
            if (val < 1 || val > 60 || c != '\0')
            {
                logmsg ("Invalid value for gcint=\n");
                return -1;
            }
            else
            {
                cckdblk.gcwait = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "gcparm") == 0)
        {
            if (val < -8 || val > 8 || c != '\0')
            {
                logmsg ("Invalid value for gcparm=\n");
                return -1;
            }
            else
            {
                cckdblk.gcparm = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "nostress") == 0)
        {
            if (val < 0 || val > 1 || c != '\0')
            {
                logmsg ("Invalid value for nostress=\n");
                return -1;
            }
            else
            {
                cckdblk.nostress = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "freepend") == 0)
        {
            if (val < -1 || val > CCKD_MAX_FREEPEND || c != '\0')
            {
                logmsg ("Invalid value for freepend=\n");
                return -1;
            }
            else
            {
                cckdblk.freepend = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "fsync") == 0)
        {
            if (val < 0 || val > 1 || c != '\0')
            {
                logmsg ("Invalid value for fsync=\n");
                return -1;
            }
            else
            {
                cckdblk.fsync = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "trace") == 0)
        {
            if (val < 0 || val > CCKD_MAX_TRACE || c != '\0')
            {
                logmsg ("Invalid value for trace=\n");
                return -1;
            }
            else
            {
                /* Disable tracing in case it's already active */
                CCKD_TRACE *p = cckdblk.itrace;
                cckdblk.itrace = NULL;
                if (p)
                {
                    SLEEP (1);
                    cckdblk.itrace = cckdblk.itracep = cckdblk.itracex = NULL;
                    cckdblk.itracen = 0;
                    cckd_free (NULL, "trace", p);
                }

                /* Get a new trace table */
                if (val > 0)
                {
                    p = cckd_calloc (NULL, "trace", val, sizeof(CCKD_TRACE));
                    if (p)
                    {
                        cckdblk.itracen = val;
                        cckdblk.itracex = p + val;
                        cckdblk.itracep = p;
                        cckdblk.itrace  = p;
                    }
                }
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "linuxnull") == 0)
        {
            if (val < 0 || val > 1 || c != '\0')
            {
                logmsg ("Invalid value for linuxnull=\n");
                return -1;
            }
            else
            {
                cckdblk.linuxnull = val;
                opts = 1;
            }
        }
        else if (strcasecmp (kw, "gcstart") == 0)
        {
            DEVBLK *dev;
            CCKDDASD_EXT *cckd;
            TID tid;
            int flag = 0;

            cckd_lock_devchain(0);
            for (dev = cckdblk.dev1st; dev; dev = cckd->devnext)
            {
                cckd = dev->cckd_ext;
                obtain_lock (&cckd->filelock);
                if (cckd->cdevhdr[cckd->sfn].free_total)
                {
                    cckd->cdevhdr[cckd->sfn].options |= (CCKD_OPENED | CCKD_ORDWR);
                    cckd_write_chdr (dev);
                    flag = 1;
                }
                release_lock (&cckd->filelock);
            }
            cckd_unlock_devchain();
            if (flag && cckdblk.gcs < cckdblk.gcmax)
                create_thread (&tid, JOINABLE, cckd_gcol, NULL, "cckd_gcol");
        }
        else
        {
            logmsg ("cckd invalid keyword: %s\n",kw);
            if (!cmd) return -1;
            cckd_command_help ();
            op = NULL;
        }
    }

    if (cmd && opts) cckd_command_opts();
    return 0;
} /* end function cckd_command */

/*-------------------------------------------------------------------*/
/* Print internal trace                                              */
/*-------------------------------------------------------------------*/
DLL_EXPORT void cckd_print_itrace()
{
CCKD_TRACE     *i, *p;                  /* Trace table pointers      */

    if (!cckdblk.itrace) return;
    logmsg (_("HHCCD900I print_itrace\n"));
    i = cckdblk.itrace;
    cckdblk.itrace = NULL;
    SLEEP (1);
    p = cckdblk.itracep;
    if (p >= cckdblk.itracex) p = i;
    do
    {
        if (p[0] != '\0')
            logmsg ("%s", (char *)p);
        if (++p >= cckdblk.itracex) p = i;
    } while (p != cckdblk.itracep);
    memset (i, 0, cckdblk.itracen * sizeof(CCKD_TRACE));
    cckdblk.itracep = i;
    cckdblk.itrace  = i;
} /* end function cckd_print_itrace */

/*-------------------------------------------------------------------*/
/* Write internal trace entry                                        */
/*-------------------------------------------------------------------*/
void cckd_trace(DEVBLK *dev, char *msg, ...)
{
va_list         vl;
struct timeval  tv;
time_t          t;
char            tbuf[64];
CCKD_TRACE     *p;
int             l;

    if (dev && (dev->ccwtrace||dev->ccwstep))
    {
        char *bfr;
        int  sz=1024,rc;
        bfr=malloc(1024);
        va_start(vl,msg);
        while(1)
        {
            rc=vsnprintf(bfr,sz,msg,vl);
            if(rc<0)
            {
                free(bfr);
                bfr=NULL;
                break;
            }
            if(rc<sz)
            {
                break;
            }
            sz+=256;
            bfr=realloc(bfr,sz);
        }
        if(bfr)
        {
            logmsg("%4.4X:%s",dev->devnum,bfr);
        }
        va_end(vl);
    }
    /* ISW FIXME : The following code has a potential */
    /*             for heap overrun (vsprintf)        */
    if (cckdblk.itrace)
    {
        gettimeofday(&tv, NULL);
        t = tv.tv_sec;
        strcpy(tbuf, ctime(&t));
        tbuf[19] = '\0';

        va_start(vl,msg);

        p = cckdblk.itracep++;
        if (p >= cckdblk.itracex)
        {
            p = cckdblk.itrace;
            cckdblk.itracep = p + 1;
        }

        if (p)
        {
            l = sprintf ((char *)p, "%s" "." "%6.6ld %4.4X:",
                      tbuf+11, (long)(tv.tv_usec), dev ? dev->devnum : 0);
            vsprintf ((char *)p + l, msg, vl);
        }
    }
} /* end function cckd_trace */

DEVHND cckddasd_device_hndinfo = {
        &ckddasd_init_handler,          /* Device Initialisation      */
        &ckddasd_execute_ccw,           /* Device CCW execute         */
        &cckddasd_close_device,         /* Device Close               */
        &ckddasd_query_device,          /* Device Query               */
        &cckddasd_start,                /* Device Start channel pgm   */
        &cckddasd_end,                  /* Device End channel pgm     */
        &cckddasd_start,                /* Device Resume channel pgm  */
        &cckddasd_end,                  /* Device Suspend channel pgm */
        &cckd_read_track,               /* Device Read                */
        &cckd_update_track,             /* Device Write               */
        &cckd_used,                     /* Device Query used          */
        NULL,                           /* Device Reserve             */
        NULL,                           /* Device Release             */
        NULL,                           /* Device Attention           */
        NULL,                           /* Immediate CCW Codes        */
        NULL,                           /* Signal Adapter Input       */
        NULL,                           /* Signal Adapter Ouput       */
        &ckddasd_hsuspend,              /* Hercules suspend           */
        &ckddasd_hresume                /* Hercules resume            */
};

DEVHND cfbadasd_device_hndinfo = {
        &fbadasd_init_handler,          /* Device Initialisation      */
        &fbadasd_execute_ccw,           /* Device CCW execute         */
        &cckddasd_close_device,         /* Device Close               */
        &fbadasd_query_device,          /* Device Query               */
        &cckddasd_start,                /* Device Start channel pgm   */
        &cckddasd_end,                  /* Device End channel pgm     */
        &cckddasd_start,                /* Device Resume channel pgm  */
        &cckddasd_end,                  /* Device Suspend channel pgm */
        &cfba_read_block,               /* Device Read                */
        &cfba_write_block,              /* Device Write               */
        &cfba_used,                     /* Device Query used          */
        NULL,                           /* Device Reserve             */
        NULL,                           /* Device Release             */
        NULL,                           /* Device Attention           */
        NULL,                           /* Immediate CCW Codes        */
        NULL,                           /* Signal Adapter Input       */
        NULL,                           /* Signal Adapter Ouput       */
        &fbadasd_hsuspend,              /* Hercules suspend           */
        &fbadasd_hresume                /* Hercules resume            */
};
