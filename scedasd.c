/* SCEDASD.C    (c) Copyright Jan Jaeger, 2001-2009                  */
/*              Service Control Element DASD I/O functions           */

#include "hstdinc.h"

#define _HENGINE_DLL_

#include "hercules.h"
#include "opcode.h"

#if !defined(_SCEDASD_C)
#define _SCEDASD_C

static TID     scedio_tid;             /* Thread id of the i/o driver
                                                                     */
static char *sce_basedir = NULL;


char *get_sce_dir()
{
    return sce_basedir;
}


void set_sce_dir(char *path)
{
char realdir[MAX_PATH];
char tempdir[MAX_PATH];

    if(sce_basedir)
    {
        free(sce_basedir);
        sce_basedir = NULL;
    }

    if(!path)
        sce_basedir = NULL;
    else
        if(!realpath(path,tempdir))
        {
            logmsg(_("HHCSC011E set_sce_dir: %s: %s\n"),path,strerror(errno));
            sce_basedir = NULL;
        }
        else
        {
            hostpath(realdir, tempdir, sizeof(realdir));
            strlcat(realdir,"/",sizeof(realdir));
            sce_basedir = strdup(realdir);
        }
}


static char *set_sce_basedir(char *path)
{
char *basedir;
char realdir[MAX_PATH];
char tempdir[MAX_PATH];

    if(sce_basedir)
    {
        free(sce_basedir);
        sce_basedir = NULL;
    }

    if(!realpath(path,tempdir))
    {
        logmsg(_("HHCSC012E set_sce_basedir: %s: %s\n"),path,strerror(errno));
        sce_basedir = NULL;
        return NULL;
    }
    hostpath(realdir, tempdir, sizeof(realdir));

    if((basedir = strrchr(realdir,'/')))
    {
        *(++basedir) = '\0';
        sce_basedir = strdup(realdir);
        return (basedir = strrchr(path,'/')) ? ++basedir : path;
    }
    else
    {
        sce_basedir = NULL; 
        return path;
    }
}


static char *check_sce_filepath(const char *path, char *fullpath)
{
char temppath[MAX_PATH];
char tempreal[MAX_PATH];

    /* Return file access error if no basedir has been set */
    if(!sce_basedir)
    {
        strlcpy(fullpath,path,sizeof(temppath));
        errno = EACCES;
        return NULL;
    }

    /* Establish the full path of the file we are trying to access */
    strlcpy(temppath,sce_basedir,sizeof(temppath));
    strlcat(temppath,path,sizeof(temppath));
    
    if(!realpath(temppath,tempreal))
    {
        hostpath(fullpath, tempreal, sizeof(temppath));
        if(strncmp( sce_basedir, fullpath, strlen(sce_basedir)))
            errno = EACCES;
        return NULL;
    }
    
    hostpath(fullpath, tempreal, sizeof(temppath));
    if(strncmp( sce_basedir, fullpath, strlen(sce_basedir)))
    {
        errno = EACCES;
        return NULL;
    }

    return fullpath;
}

#endif /* !defined(_SCEDASD_C) */


/*-------------------------------------------------------------------*/
/* function load_hmc simulates the load from the service processor   */
/*   the filename pointed to is a descriptor file which has the      */
/*   following format:                                               */
/*                                                                   */
/*   '*' in col 1 is comment                                         */
/*   core image file followed by address where it should be loaded   */
/*                                                                   */
/* For example:                                                      */
/*                                                                   */
/* * Linux/390 cdrom boot image                                      */
/* boot_images/tapeipl.ikr 0x00000000                                */
/* boot_images/initrd 0x00800000                                     */
/* boot_images/parmfile 0x00010480                                   */
/*                                                                   */
/* The location of the image files is relative to the location of    */
/* the descriptor file.                         Jan Jaeger 10-11-01  */
/*                                                                   */
/*-------------------------------------------------------------------*/
int ARCH_DEP(load_hmc) (char *fname, int cpu, int clear)
{
REGS   *regs;                           /* -> Regs                   */
FILE   *fp;
char    inputbuff[MAX_PATH];
char   *inputline;
char    filename[MAX_PATH];                 /* filename of image file    */
char    pathname[MAX_PATH];                 /* pathname of image file    */
U32     fileaddr;
int     rc = 0;                         /* Return codes (work)       */

    /* Get started */
    if (ARCH_DEP(common_load_begin) (cpu, clear) != 0)
        return -1;

    /* The actual IPL proper starts here... */

    regs = sysblk.regs[cpu];    /* Point to IPL CPU's registers */

    if(fname == NULL)                   /* Default ipl from DASD     */
        fname = "HERCULES.ins";         /*   from HERCULES.ins       */

    hostpath(pathname, fname, sizeof(pathname));

    if(!(fname = set_sce_basedir(pathname)))
        return -1;

    /* Construct and check full pathname */
    if(!check_sce_filepath(fname,filename))
    {
        logmsg(_("HHCSC001E Load from %s failed: %s\n"),fname,strerror(errno));
        return -1;
    }
    
    fp = fopen(filename, "r");
    if(fp == NULL)
    {
        logmsg(_("HHCSC002E Load from %s failed: %s\n"),fname,strerror(errno));
        return -1;
    }

    do
    {
        inputline = fgets(inputbuff,sizeof(inputbuff),fp);

#if !defined(_MSVC_)
        if(inputline && *inputline == 0x1a)
            inputline = NULL;
#endif /*!defined(_MSVC_)*/

        if(inputline)
        {
            rc = sscanf(inputline,"%" MSTRING(MAX_PATH) "s %i",filename,&fileaddr);
        }

        /* If no load address was found load to location zero */
        if(inputline && rc < 2)
            fileaddr = 0;

        if(inputline && rc > 0 && *filename != '*' && *filename != '#')
        {
            hostpath(pathname, filename, sizeof(pathname));

            /* Construct and check full pathname */
            if(!check_sce_filepath(pathname,filename))
            {
                logmsg(_("HHCSC003E Load from %s failed: %s\n"),pathname,strerror(errno));
                return -1;
            }

            if( ARCH_DEP(load_main) (filename, fileaddr) < 0 )
            {
                fclose(fp);
                HDC1(debug_cpu_state, regs);
                return -1;
            }
            sysblk.main_clear = sysblk.xpnd_clear = 0;
        }
    } while(inputline);
    fclose(fp);

    /* Finish up... */
    return ARCH_DEP(common_load_finish) (regs);

} /* end function load_hmc */


/*-------------------------------------------------------------------*/
/* Function to Load (read) specified file into absolute main storage */
/*-------------------------------------------------------------------*/
int ARCH_DEP(load_main) (char *fname, RADR startloc)
{
int fd;
int len;
int rc = 0;
RADR pageaddr;
U32  pagesize;

    fd = hopen(fname, O_RDONLY|O_BINARY);
    if (fd < 0)
    {
        if(errno != ENOENT)
            logmsg(_("HHCSC031E load_main: %s: %s\n"), fname, strerror(errno));
        return fd;
    }

    pagesize = PAGEFRAME_PAGESIZE - (startloc & PAGEFRAME_BYTEMASK);
    pageaddr = startloc;

    for( ; ; ) {
        if (pageaddr >= sysblk.mainsize)
        {
            logmsg(_("HHCSC032W load_main: terminated at end of mainstor\n"));
            close(fd);
            return rc;
        }

        len = read(fd, sysblk.mainstor + pageaddr, pagesize);
        if (len > 0)
        {
            STORAGE_KEY(pageaddr, &sysblk) |= STORKEY_REF|STORKEY_CHANGE;
            rc += len;
        }

        if (len < (int)pagesize)
        {
            close(fd);
            return rc;
        }

        pageaddr += PAGEFRAME_PAGESIZE;
        pageaddr &= PAGEFRAME_PAGEMASK;
        pagesize  = PAGEFRAME_PAGESIZE;
    }

} /* end function load_main */


#if defined(FEATURE_SCEDIO)

/*-------------------------------------------------------------------*/
/* Function to write to a file on the service processor disk         */
/*-------------------------------------------------------------------*/
static S64 ARCH_DEP(write_file)(char *fname, int mode, CREG sto, S64 size)
{
int fd, nwrite;
U64 totwrite = 0;

    fd = hopen(fname, mode |O_WRONLY|O_BINARY,
            S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
    if (fd < 0)
    {
        logmsg (_("HHCSC041E %s open error: %s\n"), fname, strerror(errno));
        return -1;
    }

#if defined(FEATURE_ESAME)
    sto &= ASCE_TO;
#else /*!defined(FEATURE_ESAME)*/
    sto &= STD_STO;
#endif /*!defined(FEATURE_ESAME)*/

    for( ; size > 0 ; sto += sizeof(sto))
    {
#if defined(FEATURE_ESAME)
    DBLWRD *ste;
#else /*!defined(FEATURE_ESAME)*/
    FWORD *ste;
#endif /*!defined(FEATURE_ESAME)*/
    CREG pto, pti;

        /* Fetch segment table entry and calc Page Table Origin */
        if( sto >= sysblk.mainsize)
            goto eof;
#if defined(FEATURE_ESAME)
        ste = (DBLWRD*)(sysblk.mainstor + sto);
#else /*!defined(FEATURE_ESAME)*/
        ste = (FWORD*)(sysblk.mainstor + sto);
#endif /*!defined(FEATURE_ESAME)*/
        FETCH_W(pto, ste);
        STORAGE_KEY(sto, &sysblk) |= (STORKEY_REF);
        if( pto & SEGTAB_INVALID )
            goto eof;
#if defined(FEATURE_ESAME)
        pto &= ZSEGTAB_PTO;
#else /*!defined(FEATURE_ESAME)*/
        pto &= SEGTAB_PTO;
#endif /*!defined(FEATURE_ESAME)*/

        for(pti = 0; pti < 256 && size > 0; pti++, pto += sizeof(pto))
        {
#if defined(FEATURE_ESAME)
        DBLWRD *pte;
#else /*!defined(FEATURE_ESAME)*/
        FWORD *pte;
#endif /*!defined(FEATURE_ESAME)*/
        CREG pgo;
        BYTE *page;

            /* Fetch Page Table Entry to get page origin */
            if( pto >= sysblk.mainsize)
                goto eof;

#if defined(FEATURE_ESAME)
            pte = (DBLWRD*)(sysblk.mainstor + pto);
#else /*!defined(FEATURE_ESAME)*/
            pte = (FWORD*)(sysblk.mainstor + pto);
#endif /*!defined(FEATURE_ESAME)*/
            FETCH_W(pgo, pte);
            STORAGE_KEY(pto, &sysblk) |= (STORKEY_REF);
            if( !(pgo & PAGETAB_INVALID) )
            {
#if defined(FEATURE_ESAME)
                pgo &= ZPGETAB_PFRA;
#else /*!defined(FEATURE_ESAME)*/
                pgo &= PAGETAB_PFRA;
#endif /*!defined(FEATURE_ESAME)*/

                /* Write page to SCE disk */
                if( pgo >= sysblk.mainsize)
                    goto eof;

                page = sysblk.mainstor + pgo;
                nwrite = write(fd, page, STORAGE_KEY_PAGESIZE);
                totwrite += nwrite;
                if( nwrite != STORAGE_KEY_PAGESIZE )
                    goto eof;

                STORAGE_KEY(pgo, &sysblk) |= (STORKEY_REF);
            }
            size -= STORAGE_KEY_PAGESIZE;
        }
    }
eof:
    close(fd);
    return totwrite;
}


/*-------------------------------------------------------------------*/
/* Function to read from a file on the service processor disk        */
/*-------------------------------------------------------------------*/
static S64 ARCH_DEP(read_file)(char *fname, CREG sto, S64 seek, S64 size)
{
int fd, nread;
U64 totread = 0;

    fd = hopen(fname, O_RDONLY|O_BINARY);
    if (fd < 0)
    {
        if(errno != ENOENT)
            logmsg (_("HHCSC051E %s open error: %s\n"), fname, strerror(errno));
        return -1;
    }

    if(lseek(fd, (off_t)seek, SEEK_SET) == (off_t)seek)
    {
#if defined(FEATURE_ESAME)
        sto &= ASCE_TO;
#else /*!defined(FEATURE_ESAME)*/
        sto &= STD_STO;
#endif /*!defined(FEATURE_ESAME)*/

        for( ; size > 0 ; sto += sizeof(sto))
        {
#if defined(FEATURE_ESAME)
        DBLWRD *ste;
#else /*!defined(FEATURE_ESAME)*/
        FWORD *ste;
#endif /*!defined(FEATURE_ESAME)*/
        CREG pto, pti;

            /* Fetch segment table entry and calc Page Table Origin */
            if( sto >= sysblk.mainsize)
                goto eof;
#if defined(FEATURE_ESAME)
            ste = (DBLWRD*)(sysblk.mainstor + sto);
#else /*!defined(FEATURE_ESAME)*/
            ste = (FWORD*)(sysblk.mainstor + sto);
#endif /*!defined(FEATURE_ESAME)*/
            FETCH_W(pto, ste);
            STORAGE_KEY(sto, &sysblk) |= (STORKEY_REF);
            if( pto & SEGTAB_INVALID )
                goto eof;
#if defined(FEATURE_ESAME)
            pto &= ZSEGTAB_PTO;
#else /*!defined(FEATURE_ESAME)*/
            pto &= SEGTAB_PTO;
#endif /*!defined(FEATURE_ESAME)*/

            for(pti = 0; pti < 256 && size > 0; pti++, pto += sizeof(pto))
            {
#if defined(FEATURE_ESAME)
            DBLWRD *pte;
#else /*!defined(FEATURE_ESAME)*/
            FWORD *pte;
#endif /*!defined(FEATURE_ESAME)*/
            CREG pgo;
            BYTE *page;

                /* Fetch Page Table Entry to get page origin */
                if( pto >= sysblk.mainsize)
                    goto eof;
#if defined(FEATURE_ESAME)
                pte = (DBLWRD*)(sysblk.mainstor + pto);
#else /*!defined(FEATURE_ESAME)*/
                pte = (FWORD*)(sysblk.mainstor + pto);
#endif /*!defined(FEATURE_ESAME)*/
                FETCH_W(pgo, pte);
                STORAGE_KEY(pto, &sysblk) |= (STORKEY_REF);
                if( pgo & PAGETAB_INVALID )
                    goto eof;
#if defined(FEATURE_ESAME)
                pgo &= ZPGETAB_PFRA;
#else /*!defined(FEATURE_ESAME)*/
                pgo &= PAGETAB_PFRA;
#endif /*!defined(FEATURE_ESAME)*/

                /* Read page into main storage */
                if( pgo >= sysblk.mainsize)
                    goto eof;
                page = sysblk.mainstor + pgo;
                nread = read(fd, page, STORAGE_KEY_PAGESIZE);
                totread += nread;
                size -= nread;
                if( nread != STORAGE_KEY_PAGESIZE )
                    goto eof;
                STORAGE_KEY(pgo, &sysblk) |= (STORKEY_REF|STORKEY_CHANGE);
            }
        }
    }
eof:
    close(fd);
    return totread;
}


static int ARCH_DEP(scedio_ior)(SCCB_SCEDIOR_BK *scedior_bk)
{
U32  origin;
char image[9];
S32  size;
unsigned int i;
char filename[MAX_PATH];

    FETCH_FW(origin,scedior_bk->origin);

    /* Convert image filename to null terminated ascii string */
    for(i = 0; i < sizeof(image)-1 && scedior_bk->image[i] != 0x40; i++)
        image[i] = guest_to_host((int)scedior_bk->image[i]);
    image[i] = '\0';

    /* Ensure file access is allowed and within specified directory */
    if(!check_sce_filepath(image,filename))
    {
        if(errno != ENOENT)
            logmsg (_("HHCSC101E access error: %s: %s\n"), image, strerror(errno));
        return FALSE;
    }

    size = ARCH_DEP(load_main)(filename,origin);

    return (size >= 0) ? TRUE : FALSE;
}


static int ARCH_DEP(scedio_iov)(SCCB_SCEDIOV_BK *scediov_bk)
{
S64     seek;
S64     length;
S64     totread, totwrite;
U64     sto;
char    fname[MAX_PATH];

    switch(scediov_bk->type) {

    case SCCB_SCEDIOV_TYPE_INIT:
        return TRUE;
        break;


    case SCCB_SCEDIOV_TYPE_READ:
        /* Ensure file access is allowed and within specified directory */
        if(!check_sce_filepath((char*)scediov_bk->filename,fname))
        {
            if(errno != ENOENT)
                logmsg (_("HHCSC201E access error: %s: %s\n"), fname, strerror(errno));
            return FALSE;
        }
        FETCH_DW(sto,scediov_bk->sto);
        FETCH_DW(seek,scediov_bk->seek);
        FETCH_DW(length,scediov_bk->length);

        totread = ARCH_DEP(read_file)(fname, sto, seek, length);

        if(totread > 0)
        {
            STORE_DW(scediov_bk->length,totread);

            if(totread == length)
                STORE_DW(scediov_bk->ncomp,0);
            else
                STORE_DW(scediov_bk->ncomp,seek+totread);

            return TRUE;
        }
        else
            return FALSE;

        break;


    case SCCB_SCEDIOV_TYPE_CREATE:
    case SCCB_SCEDIOV_TYPE_APPEND:
        /* Ensure file access is allowed and within specified directory */
        if(!check_sce_filepath((char*)scediov_bk->filename,fname))
        {
            if(errno != ENOENT)
                logmsg (_("HHCSC202I access error: %s: %s\n"), fname, strerror(errno));

            /* A file not found error may be expected for a create request */
            if(!(errno == ENOENT && scediov_bk->type == SCCB_SCEDIOV_TYPE_CREATE))
                return FALSE;
        }
        FETCH_DW(sto,scediov_bk->sto);
        FETCH_DW(seek,scediov_bk->seek);
        FETCH_DW(length,scediov_bk->length);

        totwrite = ARCH_DEP(write_file)(fname,
          ((scediov_bk->type == SCCB_SCEDIOV_TYPE_CREATE) ? (O_CREAT|O_TRUNC) : O_APPEND), sto, length);

        if(totwrite >= 0)
        {
            STORE_DW(scediov_bk->ncomp,totwrite);
            return TRUE;
        }
        else
            return FALSE;
        break;


    default:
        PTT(PTT_CL_ERR,"*SERVC",(U32)scediov_bk->type,(U32)scediov_bk->flag1,scediov_bk->flag2);
        return FALSE;

    }
}


/*-------------------------------------------------------------------*/
/* Thread to perform service processor I/O functions                 */
/*-------------------------------------------------------------------*/
static void *ARCH_DEP(scedio_thread)(SCCB_SCEDIO_BK *scedio_bk)
{
SCCB_SCEDIOV_BK *scediov_bk;
SCCB_SCEDIOR_BK *scedior_bk;

    switch(scedio_bk->flag1) {

    case SCCB_SCEDIO_FLG1_IOV:
        scediov_bk = (SCCB_SCEDIOV_BK*)(scedio_bk + 1);
        if( ARCH_DEP(scedio_iov)(scediov_bk) )
            scedio_bk->flag3 |= SCCB_SCEDIO_FLG3_COMPLETE;
        else
            scedio_bk->flag3 &= ~SCCB_SCEDIO_FLG3_COMPLETE;
        break;

    case SCCB_SCEDIO_FLG1_IOR:
        scedior_bk = (SCCB_SCEDIOR_BK*)(scedio_bk + 1);
        if( ARCH_DEP(scedio_ior)(scedior_bk) )
            scedio_bk->flag3 |= SCCB_SCEDIO_FLG3_COMPLETE;
        else
            scedio_bk->flag3 &= ~SCCB_SCEDIO_FLG3_COMPLETE;
        break;

    default:
        PTT(PTT_CL_ERR,"*SERVC",(U32)scedio_bk->flag0,(U32)scedio_bk->flag1,scedio_bk->flag3);
    }


    OBTAIN_INTLOCK(NULL);

    while(IS_IC_SERVSIG)
    {
        RELEASE_INTLOCK(NULL);
        sched_yield();
        OBTAIN_INTLOCK(NULL);
    }

    sclp_attention(SCCB_EVD_TYPE_SCEDIO);

    scedio_tid = 0;

    RELEASE_INTLOCK(NULL);
    return NULL;
}


/*-------------------------------------------------------------------*/
/* Function to interface with the service processor I/O thread       */
/*-------------------------------------------------------------------*/
static int ARCH_DEP(scedio_request)(U32 sclp_command, SCCB_EVD_HDR *evd_hdr)
{
SCCB_SCEDIO_BK  *scedio_bk = (SCCB_SCEDIO_BK*)(evd_hdr + 1);
SCCB_SCEDIOV_BK *scediov_bk;
SCCB_SCEDIOR_BK *scedior_bk;


static struct {
    SCCB_SCEDIO_BK scedio_bk;
    union {
        SCCB_SCEDIOV_BK v;
        SCCB_SCEDIOR_BK r;
        } io;
    } static_scedio_bk ;

static int scedio_pending;

    if(sclp_command == SCLP_READ_EVENT_DATA)
    {
    int pending_req = scedio_pending;
    U16 evd_len;

        /* Return no data if the scedio thread is still active */
        if(scedio_tid)
            return 0;

        /* Update the scedio_bk copy in the SCCB */
        if(scedio_pending)
        {
            /* Zero all fields */
            memset (evd_hdr, 0, sizeof(SCCB_EVD_HDR));

            /* Set type in event header */
            evd_hdr->type = SCCB_EVD_TYPE_SCEDIO;

            /* Store scedio header */
            *scedio_bk = static_scedio_bk.scedio_bk;

            /* Calculate event response length */
            evd_len = sizeof(SCCB_EVD_HDR) + sizeof(SCCB_SCEDIO_BK);

            switch(scedio_bk->flag1) {
            case SCCB_SCEDIO_FLG1_IOR:
                scedior_bk = (SCCB_SCEDIOR_BK*)(scedio_bk + 1);
                *scedior_bk = static_scedio_bk.io.r;
                evd_len += sizeof(SCCB_SCEDIOR_BK);
                break;
            case SCCB_SCEDIO_FLG1_IOV:
                scediov_bk = (SCCB_SCEDIOV_BK*)(scedio_bk + 1);
                *scediov_bk = static_scedio_bk.io.v ;
                evd_len += sizeof(SCCB_SCEDIOV_BK);
                break;
            default:
                PTT(PTT_CL_ERR,"*SERVC",(U32)evd_hdr->type,(U32)scedio_bk->flag1,scedio_bk->flag3);
            }

            /* Set length in event header */
            STORE_HW(evd_hdr->totlen, evd_len);
        }


        /* Reset the pending flag */
        scedio_pending = 0;

        /* Return true if a request was pending */
        return pending_req;

    }
    else
    {
#if !defined(NO_SIGABEND_HANDLER)
        switch(scedio_bk->flag1) {
            case SCCB_SCEDIO_FLG1_IOV:
                scediov_bk = (SCCB_SCEDIOV_BK*)(scedio_bk + 1);
                switch(scediov_bk->type) {
                    case SCCB_SCEDIOV_TYPE_INIT:
                        /* Kill the scedio thread if it is active */
                        if( scedio_tid )
                        {
                            OBTAIN_INTLOCK(NULL);
                            signal_thread(scedio_tid, SIGKILL);
                            scedio_tid = 0;
                            scedio_pending = 0;
                            RELEASE_INTLOCK(NULL);
                        }
                        break;
                }
                break;
        }
#endif

        /* Take a copy of the scedio_bk in the SCCB */
        static_scedio_bk.scedio_bk = *scedio_bk;
        switch(scedio_bk->flag1) {
        case SCCB_SCEDIO_FLG1_IOR:
            scedior_bk = (SCCB_SCEDIOR_BK*)(scedio_bk + 1);
            static_scedio_bk.io.r = *scedior_bk;
            break;
        case SCCB_SCEDIO_FLG1_IOV:
            scediov_bk = (SCCB_SCEDIOV_BK*)(scedio_bk + 1);
            static_scedio_bk.io.v = *scediov_bk;
            break;
        default:
            PTT(PTT_CL_ERR,"*SERVC",(U32)evd_hdr->type,(U32)scedio_bk->flag1,scedio_bk->flag3);
        }

        /* Create the scedio thread */
        if( create_thread(&scedio_tid, &sysblk.detattr,
            ARCH_DEP(scedio_thread), &static_scedio_bk, "scedio_thread") )
            return -1;

        scedio_pending = 1;

    }

    return 0;
}


/*-------------------------------------------------------------------*/
/* Function to request service processor I/O                         */
/*-------------------------------------------------------------------*/
void ARCH_DEP(sclp_scedio_request) (SCCB_HEADER *sccb)
{
SCCB_EVD_HDR    *evd_hdr = (SCCB_EVD_HDR*)(sccb + 1);

    if( ARCH_DEP(scedio_request)(SCLP_WRITE_EVENT_DATA, evd_hdr) )
    {
        /* Set response code X'0040' in SCCB header */
        sccb->reas = SCCB_REAS_NONE;
        sccb->resp = SCCB_RESP_BACKOUT;
    }
    else
    {
        /* Set response code X'0020' in SCCB header */
        sccb->reas = SCCB_REAS_NONE;
        sccb->resp = SCCB_RESP_COMPLETE;
    }

    /* Indicate Event Processed */
    evd_hdr->flag |= SCCB_EVD_FLAG_PROC;

}


/*-------------------------------------------------------------------*/
/* Function to read service processor I/O event data                 */
/*-------------------------------------------------------------------*/
void ARCH_DEP(sclp_scedio_event) (SCCB_HEADER *sccb)
{
SCCB_EVD_HDR    *evd_hdr = (SCCB_EVD_HDR*)(sccb + 1);
U16 sccb_len;
U16 evd_len;

    if( ARCH_DEP(scedio_request)(SCLP_READ_EVENT_DATA, evd_hdr) )
    {
        /* Update SCCB length field if variable request */
        if (sccb->type & SCCB_TYPE_VARIABLE)
        {
            FETCH_HW(evd_len, evd_hdr->totlen);
            sccb_len = evd_len + sizeof(SCCB_HEADER);
            STORE_HW(sccb->length, sccb_len);
            sccb->type &= ~SCCB_TYPE_VARIABLE;
        }

        /* Set response code X'0020' in SCCB header */
        sccb->reas = SCCB_REAS_NONE;
        sccb->resp = SCCB_RESP_COMPLETE;
    }

}

#endif /*defined(FEATURE_SCEDIO)*/


#if !defined(_GEN_ARCH)

#if defined(_ARCHMODE2)
 #define  _GEN_ARCH _ARCHMODE2
 #include "scedasd.c"
#endif

#if defined(_ARCHMODE3)
 #undef   _GEN_ARCH
 #define  _GEN_ARCH _ARCHMODE3
 #include "scedasd.c"
#endif


/*-------------------------------------------------------------------*/
/*  Service Processor Load    (load/ipl from the specified file)     */
/*-------------------------------------------------------------------*/

int load_hmc (char *fname, int cpu, int clear)
{
    switch(sysblk.arch_mode) {
#if defined(_370)
        case ARCH_370:
            return s370_load_hmc (fname, cpu, clear);
#endif
#if defined(_390)
        case ARCH_390:
            return s390_load_hmc (fname, cpu, clear);
#endif
#if defined(_900)
        case ARCH_900:
            /* z/Arch always starts out in ESA390 mode */
            return s390_load_hmc (fname, cpu, clear);
#endif
    }
    return -1;
}


/*-------------------------------------------------------------------*/
/* Load/Read specified file into absolute main storage               */
/*-------------------------------------------------------------------*/
int load_main (char *fname, RADR startloc)
{
    switch(sysblk.arch_mode) {
#if defined(_370)
        case ARCH_370:
            return s370_load_main (fname, startloc);
#endif
#if defined(_390)
        case ARCH_390:
            return s390_load_main (fname, startloc);
#endif
#if defined(_900)
        case ARCH_900:
            return z900_load_main (fname, startloc);
#endif
    }
    return -1;
}

#endif /*!defined(_GEN_ARCH)*/
