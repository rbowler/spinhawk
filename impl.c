/* IMPL.C       (c) Copyright Roger Bowler, 1999-2009                */
/*              Hercules Initialization Module                       */

/*-------------------------------------------------------------------*/
/* This module initializes the Hercules S/370 or ESA/390 emulator.   */
/* It builds the system configuration blocks, creates threads for    */
/* central processors, HTTP server, logger task and activates the    */
/* control panel which runs under the main thread when in foreground */
/* mode.                                                             */
/*-------------------------------------------------------------------*/

#include "hstdinc.h"

#define _IMPL_C_
#define _HENGINE_DLL_

#include "hercules.h"
#include "opcode.h"
#include "devtype.h"
#include "hconsole.h"
#include "herc_getopt.h"
#include "hostinfo.h"
#include "history.h"

/* (delayed_exit function defined in config.c) */
extern void delayed_exit (int exit_code);

/* forward define process_script_file (ISW20030220-3) */
int process_script_file(char *,int);

static LOGCALLBACK  log_callback=NULL;

/*-------------------------------------------------------------------*/
/* Register a LOG callback                                           */
/*-------------------------------------------------------------------*/
DLL_EXPORT void registerLogCallback(LOGCALLBACK lcb)
{
    log_callback=lcb;
}

/*-------------------------------------------------------------------*/
/* Signal handler for SIGINT signal                                  */
/*-------------------------------------------------------------------*/
static void sigint_handler (int signo)
{
//  logmsg ("impl.c: sigint handler entered for thread %lu\n",/*debug*/
//          thread_id());                                     /*debug*/

    UNREFERENCED(signo);

    signal(SIGINT, sigint_handler);
    /* Ignore signal unless presented on console thread */
    if ( !equal_threads( thread_id(), sysblk.cnsltid ) )
        return;

    /* Exit if previous SIGINT request was not actioned */
    if (sysblk.sigintreq)
    {
        /* Release the configuration */
        release_config();
        delayed_exit(1);
    }

    /* Set SIGINT request pending flag */
    sysblk.sigintreq = 1;

    /* Activate instruction stepping */
    sysblk.inststep = 1;
    SET_IC_TRACE;
    return;
} /* end function sigint_handler */

/*-------------------------------------------------------------------*/
/* Signal handler for SIGTERM signal                                 */
/*-------------------------------------------------------------------*/
static void sigterm_handler (int signo)
{
//  logmsg ("impl.c: sigterm handler entered for thread %lu\n",/*debug*/
//          thread_id());                                      /*debug*/

    UNREFERENCED(signo);

    signal(SIGTERM, sigterm_handler);
    /* Ignore signal unless presented on main program (impl) thread */
    if ( !equal_threads( thread_id(), sysblk.impltid ) )
        return;

    /* Initiate system shutdown */
    do_shutdown();

    return;
} /* end function sigterm_handler */

/*-------------------------------------------------------------------*/
/* Signal handler for SIGWINCH signal                                 */
/*-------------------------------------------------------------------*/
static void sigwinch_handler (int signo)
{
//  logmsg ("impl.c: sigwinch handler entered for thread %lu\n",/*debug*/
//          thread_id());                                       /*debug*/

    UNREFERENCED(signo);

     window_changed = 1;

    return;
} /* end function sigwinch_handler */

#if defined( _MSVC_ )

/*-------------------------------------------------------------------*/
/* Signal handler for Windows signals                                */
/*-------------------------------------------------------------------*/
BOOL WINAPI console_ctrl_handler (DWORD signo)
{
    int i;

    SetConsoleCtrlHandler(console_ctrl_handler, FALSE);   // turn handler off while processing

    switch ( signo )
    {
        case CTRL_BREAK_EVENT:
            logmsg(_("HHCIN050ICtrl-Break intercepted. Interrupt Key depressed simulated.\n"));

            OBTAIN_INTLOCK(NULL);

            ON_IC_INTKEY;

            /* Signal waiting CPUs that an interrupt is pending */
            WAKEUP_CPUS_MASK (sysblk.waiting_mask);

            RELEASE_INTLOCK(NULL);

            SetConsoleCtrlHandler(console_ctrl_handler, TRUE);  // reset handler
            return TRUE;
            break;
        case CTRL_C_EVENT:
            logmsg(_("HHCIN022I Ctrl-C intercepted\n"));
            SetConsoleCtrlHandler(console_ctrl_handler, TRUE);  // reset handler
            return TRUE;
            break;
        case CTRL_CLOSE_EVENT:
        case CTRL_SHUTDOWN_EVENT:
        case CTRL_LOGOFF_EVENT:
            if ( !sysblk.shutdown )  // (system shutdown not initiated)
            {
                logmsg(_("HHCIN021ICLOSE Event received, SHUTDOWN Immediate starting...\n"));
                sysblk.shutimmed = TRUE;
                do_shutdown();

//                logmsg("%s(%d): return from shutdown\n", __FILE__, __LINE__ ); /* debug */

                for ( i = 0; i < 120; i++ )
                {
                    if ( sysblk.shutdown && sysblk.shutfini )
                    {
//                        logmsg("%s(%d): %d shutdown completed\n",  /* debug */
//                                __FILE__, __LINE__, i );           /* debug */
                        sleep(1);
                        break;
                    }
                    else
                    {
//                        logmsg("%s(%d): %d waiting for shutdown to complete\n",   /* debug */
//                                __FILE__, __LINE__, i );                          /* debug */
                        sleep(1);
                    }
                }

                socket_deinit();
            }
            else
            {
                logmsg(_("HHCIN023W CLOSE Event received, SHUTDOWN previously requested...\n"));
            }
            usleep(10000);
            return FALSE;
            break;
        default:
            return FALSE;
    }

} /* end function console_ctrl_handler */
#endif

#if !defined(NO_SIGABEND_HANDLER)
static void *watchdog_thread(void *arg)
{
S64 savecount[MAX_CPU_ENGINES];
int i;

    UNREFERENCED(arg);

    /* Set watchdog priority just below cpu priority
       such that it will not invalidly detect an
       inoperable cpu */
    if(sysblk.cpuprio >= 0)
        setpriority(PRIO_PROCESS, 0, sysblk.cpuprio+1);

    for (i = 0; i < MAX_CPU_ENGINES; i ++) savecount[i] = -1;

    while(!sysblk.shutdown)
    {
        for (i = 0; i < MAX_CPU; i++)
        {
//          obtain_lock (&sysblk.cpulock[i]);
            if (IS_CPU_ONLINE(i)
             && sysblk.regs[i]->cpustate == CPUSTATE_STARTED
             && (!WAITSTATE(&sysblk.regs[i]->psw)
#if defined(_FEATURE_WAITSTATE_ASSIST)
             && !(sysblk.regs[i]->sie_active && WAITSTATE(&sysblk.regs[i]->guestregs->psw))
#endif
                                           ))
            {
                /* If the cpu is running but not executing
                   instructions then it must be malfunctioning */
                if((INSTCOUNT(sysblk.regs[i]) == (U64)savecount[i])
                  && !HDC1(debug_watchdog_signal, sysblk.regs[i]) )
                {
                    /* Send signal to looping CPU */
                    signal_thread(sysblk.cputid[i], SIGUSR1);
                    savecount[i] = -1;
                }
                else
                    /* Save current instcount */
                    savecount[i] = INSTCOUNT(sysblk.regs[i]);
            }
            else
                /* mark savecount invalid as CPU not in running state */
                savecount[i] = -1;
//          release_lock (&sysblk.cpulock[i]);
        }
        /* Sleep for 20 seconds */
        SLEEP(20);
    }

    return NULL;
}
#endif /*!defined(NO_SIGABEND_HANDLER)*/

void *log_do_callback(void *dummy)
{
    char *msgbuf;
    int msgcnt = -1,msgnum;
    UNREFERENCED(dummy);
    while(msgcnt)
    {
        if((msgcnt = log_read(&msgbuf, &msgnum, LOG_BLOCK)))
        {
            log_callback(msgbuf,msgcnt);
        }
    }
    return(NULL);
}

DLL_EXPORT  COMMANDHANDLER getCommandHandler(void)
{
    return(panel_command);
}

/*-------------------------------------------------------------------*/
/* Process .RC file thread                                           */
/*-------------------------------------------------------------------*/

void* process_rc_file (void* dummy)
{
char   *rcname;                         /* hercules.rc name pointer  */
int     is_default_rc  = 0;             /* 1 == default name used    */
int     numcpu         = 0;             /* #of ONLINE & STOPPED CPUs */
int     i;                              /* (work)                    */

    UNREFERENCED(dummy);

    /* Wait for all installed/configured CPUs to
       come ONLINE and enter the STOPPED state */

    OBTAIN_INTLOCK(NULL);

    for (;;)
    {
        numcpu = 0;
        for (i = 0; i < MAX_CPU_ENGINES; i++)
            if (IS_CPU_ONLINE(i) &&
                CPUSTATE_STOPPED == sysblk.regs[i]->cpustate)
                numcpu++;
        if (numcpu == sysblk.numcpu)
            break;
        RELEASE_INTLOCK(NULL);
        usleep( 10 * 1000 );
        OBTAIN_INTLOCK(NULL);
    }

    RELEASE_INTLOCK(NULL);

    /* Wait for panel thread to engage */

    while (!sysblk.panel_init)
        usleep( 10 * 1000 );

    /* Obtain the name of the hercules.rc file or default */

    if (!(rcname = getenv("HERCULES_RC")))
    {
        rcname = "hercules.rc";
        is_default_rc = 1;
    }

#if defined(OPTION_HAO)
    /* Initialize the Hercules Automatic Operator */

    if ( !hao_initialize() )
        logmsg(_("HHCIN004S Cannot create HAO thread: %s\n"),
                strerror(errno));
#endif /* defined(OPTION_HAO) */

    /* Run the script processor for this file */

    if (process_script_file(rcname,1) != 0)
        if (ENOENT == errno)
            if (!is_default_rc)
                logmsg(_("HHCPN995E .RC file \"%s\" not found.\n"),
                    rcname);
        // (else error message already issued)

    return NULL;
}


/*-------------------------------------------------------------------*/
/* IMPL main entry point                                             */
/*-------------------------------------------------------------------*/
DLL_EXPORT int impl(int argc, char *argv[])
{
char   *cfgfile;                        /* -> Configuration filename */
int     c;                              /* Work area for getopt      */
int     arg_error = 0;                  /* 1=Invalid arguments       */
char   *msgbuf;                         /*                           */
int     msgnum;                         /*                           */
int     msgcnt;                         /*                           */
TID     rctid;                          /* RC file thread identifier */
TID     logcbtid;                       /* RC file thread identifier */

    SET_THREAD_NAME("impl");

    /* Initialize 'hostinfo' BEFORE display_version is called */
    init_hostinfo( &hostinfo );

#ifdef _MSVC_
    /* Initialize sockets package */
    VERIFY( socket_init() == 0 );
#endif

    /* Ensure hdl_shut is called in case of shutdown
       hdl_shut will ensure entries are only called once */
    atexit(hdl_shut);

    set_codepage(NULL);

    /* Clear the system configuration block */
    memset (&sysblk, 0, sizeof(SYSBLK));

    /* Save thread ID of main program */
    sysblk.impltid = thread_id();

    /* Save TOD of when we were first IMPL'ed */
    time( &sysblk.impltime );

#ifdef OPTION_MSGHLD
    /* Set the default timeout value */
    sysblk.keep_timeout_secs = 120;
#endif

    /* Initialize thread creation attributes so all of hercules
       can use them at any time when they need to create_thread
    */
    initialize_detach_attr (DETACHED);
    initialize_join_attr   (JOINABLE);

    /* Copy length for regs */
    sysblk.regs_copy_len = (int)((uintptr_t)&sysblk.dummyregs.regs_copy_end
                               - (uintptr_t)&sysblk.dummyregs);

    /* Set the daemon_mode flag indicating whether we running in
       background/daemon mode or not (meaning both stdout/stderr
       are redirected to a non-tty device). Note that this flag
       needs to be set before logger_init gets called since the
       logger_logfile_write function relies on its setting.
    */
    sysblk.daemon_mode = !isatty(STDERR_FILENO) && !isatty(STDOUT_FILENO);

    /* Initialize the logmsg pipe and associated logger thread.
       This causes all subsequent logmsg's to be redirected to
       the logger facility for handling by virtue of stdout/stderr
       being redirected to the logger facility.
    */
    logger_init();

    /* Now display the version information again after logger_init
       has been called so that either the panel display thread or the
       external gui can see the version which was previously possibly
       only displayed to the actual physical screen the first time we
       did it further above (depending on whether we're running in
       daemon_mode (external gui mode) or not). This it the call that
       the panel thread or the one the external gui actually "sees".
       The first call further above wasn't seen by either since it
       was issued before logger_init was called and thus got written
       directly to the physical screen whereas this one will be inter-
       cepted and handled by the logger facility thereby allowing the
       panel thread or external gui to "see" it and thus display it.
    */
    display_version (stdout, "Hercules ", TRUE);

#if defined(OPTION_DYNAMIC_LOAD)
    /* Initialize the hercules dynamic loader */
    hdl_main();
#endif /* defined(OPTION_DYNAMIC_LOAD) */

#ifdef EXTERNALGUI
    /* Set GUI flag if specified as final argument */
    if (argc >= 1 && strncmp(argv[argc-1],"EXTERNALGUI",11) == 0)
    {
#if defined(OPTION_DYNAMIC_LOAD)
        if (hdl_load("dyngui",HDL_LOAD_DEFAULT) != 0)
        {
            usleep(10000); /* (give logger thread time to issue
                               preceding HHCHD007E message) */
            logmsg(_("HHCIN008S DYNGUI.DLL load failed; Hercules terminated.\n"));
            delayed_exit(1);
        }
#endif /* defined(OPTION_DYNAMIC_LOAD) */
        argc--;
    }
#endif /*EXTERNALGUI*/

#if !defined(WIN32) && !defined(HAVE_STRERROR_R)
    strerror_r_init();
#endif

#if defined(OPTION_SCSI_TAPE)
    initialize_lock (&sysblk.stape_lock);
    initialize_condition (&sysblk.stape_getstat_cond);
    InitializeListHead (&sysblk.stape_mount_link);
    InitializeListHead (&sysblk.stape_status_link);
#endif /* defined(OPTION_SCSI_TAPE) */

    /* Get name of configuration file or default to hercules.cnf */
    if(!(cfgfile = getenv("HERCULES_CNF")))
        cfgfile = "hercules.cnf";

    /* Process the command line options */
    while ((c = getopt(argc, argv, "f:p:l:db:")) != EOF)
    {

        switch (c) {
        case 'f':
            cfgfile = optarg;
            break;
#if defined(OPTION_DYNAMIC_LOAD)
        case 'p':
            if(optarg)
                hdl_setpath(strdup(optarg));
            break;
        case 'l':
            {
            char *dllname, *strtok_str;
                for(dllname = strtok_r(optarg,", ",&strtok_str);
                    dllname;
                    dllname = strtok_r(NULL,", ",&strtok_str))
                    hdl_load(dllname, HDL_LOAD_DEFAULT);
            }
            break;
#endif /* defined(OPTION_DYNAMIC_LOAD) */
        case 'b':
            sysblk.logofile=optarg;
            break;
        case 'd':
            sysblk.daemon_mode = 1;
            break;
        default:
            arg_error = 1;

        } /* end switch(c) */
    } /* end while */

    if (optind < argc)
        arg_error = 1;

    /* Terminate if invalid arguments were detected */
    if (arg_error)
    {
        logmsg("usage: %s [-f config-filename] [-d] [-b logo-filename]"
#if defined(OPTION_DYNAMIC_LOAD)
                " [-p dyn-load-dir] [[-l dynmod-to-load]...]"
#endif /* defined(OPTION_DYNAMIC_LOAD) */
                " [> logfile]\n",
                argv[0]);
        delayed_exit(1);
    }

    /* Register the SIGINT handler */
    if ( signal (SIGINT, sigint_handler) == SIG_ERR )
    {
        logmsg(_("HHCIN001S Cannot register SIGINT handler: %s\n"),
                strerror(errno));
        delayed_exit(1);
    }

    /* Register the SIGTERM handler */
    if ( signal (SIGTERM, sigterm_handler) == SIG_ERR )
    {
        logmsg(_("HHCIN009S Cannot register SIGTERM handler: %s\n"),
                strerror(errno));
        delayed_exit(1);
    }

    /* Register the SIGWINCH handler */
    if ( signal (SIGWINCH, sigwinch_handler) == SIG_ERR )
    {
        logmsg(_("HHCIN010S Cannot register SIGWINCH handler: %s\n"),
                strerror(errno));
        delayed_exit(1);
    }

#if defined( _MSVC_ )
    /* Register the Window console ctrl handlers */
    if (SetConsoleCtrlHandler(console_ctrl_handler, TRUE) == FALSE)
    {
        logmsg(_("HHCIN010S Cannot register ConsoleCtrl handler: %s\n"),
                strerror(errno));
        delayed_exit(1);
    }
#endif

#if defined(HAVE_DECL_SIGPIPE) && HAVE_DECL_SIGPIPE
    /* Ignore the SIGPIPE signal, otherwise Hercules may terminate with
       Broken Pipe error if the printer driver writes to a closed pipe */
    if ( signal (SIGPIPE, SIG_IGN) == SIG_ERR )
    {
        logmsg(_("HHCIN002E Cannot suppress SIGPIPE signal: %s\n"),
                strerror(errno));
    }
#endif

#if defined( OPTION_WAKEUP_SELECT_VIA_PIPE )
    {
        int fds[2];
        initialize_lock(&sysblk.cnslpipe_lock);
        initialize_lock(&sysblk.sockpipe_lock);
        sysblk.cnslpipe_flag=0;
        sysblk.sockpipe_flag=0;
        VERIFY( create_pipe(fds) >= 0 );
        sysblk.cnslwpipe=fds[1];
        sysblk.cnslrpipe=fds[0];
        VERIFY( create_pipe(fds) >= 0 );
        sysblk.sockwpipe=fds[1];
        sysblk.sockrpipe=fds[0];
    }
#endif // defined( OPTION_WAKEUP_SELECT_VIA_PIPE )

#if !defined(NO_SIGABEND_HANDLER)
    {
    struct sigaction sa;
        sa.sa_sigaction = (void*)&sigabend_handler;
#ifdef SA_NODEFER
        sa.sa_flags = SA_NODEFER;
#else
        sa.sa_flags = 0;
#endif

        /* Explictily initialize sa_mask to its default.        @PJJ */
        sigemptyset(&sa.sa_mask);

        if( sigaction(SIGILL, &sa, NULL)
         || sigaction(SIGFPE, &sa, NULL)
         || sigaction(SIGSEGV, &sa, NULL)
         || sigaction(SIGBUS, &sa, NULL)
         || sigaction(SIGUSR1, &sa, NULL)
         || sigaction(SIGUSR2, &sa, NULL) )
        {
            logmsg(_("HHCIN003S Cannot register SIGILL/FPE/SEGV/BUS/USR "
                    "handler: %s\n"),
                    strerror(errno));
            delayed_exit(1);
        }
    }
#endif /*!defined(NO_SIGABEND_HANDLER)*/

    /* Build system configuration */
    build_config (cfgfile);

    /* System initialisation time */
    sysblk.todstart = hw_clock() << 8;

#ifdef OPTION_MIPS_COUNTING
    /* Initialize "maxrates" command reporting intervals */
    curr_int_start_time = time( NULL );
    prev_int_start_time = curr_int_start_time;
#endif

#if !defined(NO_SIGABEND_HANDLER)
    /* Start the watchdog */
    if ( create_thread (&sysblk.wdtid, DETACHED,
                        watchdog_thread, NULL, "watchdog_thread") )
    {
        logmsg(_("HHCIN004S Cannot create watchdog thread: %s\n"),
                strerror(errno));
        delayed_exit(1);
    }
#endif /*!defined(NO_SIGABEND_HANDLER)*/

#ifdef OPTION_SHARED_DEVICES
    /* Start the shared server */
    if (sysblk.shrdport)
        if ( create_thread (&sysblk.shrdtid, DETACHED,
                            shared_server, NULL, "shared_server") )
        {
            logmsg(_("HHCIN006S Cannot create shared_server thread: %s\n"),
                    strerror(errno));
            delayed_exit(1);
        }

    /* Retry pending connections */
    {
        DEVBLK *dev;
        TID     tid;

        for (dev = sysblk.firstdev; dev != NULL; dev = dev->nextdev)
            if (dev->connecting)
                if ( create_thread (&tid, DETACHED,
                           *dev->hnd->init, dev, "device connecting thread")
                   )
                {
                    logmsg(_("HHCIN007S Cannot create %4.4X connection thread: %s\n"),
                            dev->devnum, strerror(errno));
                    delayed_exit(1);
                }
    }
#endif

    /* Start up the RC file processing thread */
    create_thread(&rctid,DETACHED,
                  process_rc_file,NULL,"process_rc_file");

    if(log_callback)
    {
        // 'herclin' called us. IT'S in charge. Create its requested
        // logmsg intercept callback function and return back to it.
        create_thread(&logcbtid,DETACHED,
                      log_do_callback,NULL,"log_do_callback");
        return(0);
    }

    //---------------------------------------------------------------
    // The below functions will not return until Hercules is shutdown
    //---------------------------------------------------------------

    /* Activate the control panel */
    if(!sysblk.daemon_mode)
        panel_display ();
    else
    {
#if defined(OPTION_DYNAMIC_LOAD)
        if(daemon_task)
            daemon_task ();
        else
#endif /* defined(OPTION_DYNAMIC_LOAD) */
        {
            /* Tell RC file and HAO threads they may now proceed */
            sysblk.panel_init = 1;

            /* Retrieve messages from logger and write to stderr */
            while (1)
                if((msgcnt = log_read(&msgbuf, &msgnum, LOG_BLOCK)))
                    if(isatty(STDERR_FILENO))
                        fwrite(msgbuf,msgcnt,1,stderr);
        }
    }

    //  -----------------------------------------------------
    //      *** Hercules has been shutdown (PAST tense) ***
    //  -----------------------------------------------------

    ASSERT( sysblk.shutdown );  // (why else would we be here?!)

#ifdef _MSVC_
    SetConsoleCtrlHandler(console_ctrl_handler, FALSE);
    socket_deinit();
#endif
#ifdef DEBUG
    fprintf(stdout, _("IMPL EXIT\n"));
#endif
    fprintf(stdout, _("HHCIN099I Hercules terminated\n"));
    fflush(stdout);
    usleep(10000);
    return 0;
} /* end function main */


/*-------------------------------------------------------------------*/
/* System cleanup                                                    */
/*-------------------------------------------------------------------*/
DLL_EXPORT void system_cleanup (void)
{
//  logmsg("HHCIN950I Begin system cleanup\n");
    /*
        Currently only called by hdlmain,c's HDL_FINAL_SECTION
        after the main 'hercules' module has been unloaded, but
        that could change at some time in the future.

        The above and below logmsg's are commented out since this
        function currently doesn't do anything yet. Once it DOES
        something, they should be uncommented.
    */
//  logmsg("HHCIN959I System cleanup complete\n");
}
