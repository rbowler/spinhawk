;----------------------------------------------------------------------------
;    MODULE NAME:  Hercules-X64.MM
;
;        $Author:   USER "Dave Wade - G4UGM"  $
;      $Revision:   1.00  $
;          $Date:   14th December 2012  $
;
; DESCRIPTION
; ~~~~~~~~~~~
; This is a first attempt to build an MSI for Hercules using MAkeMsi
;----------------------------------------------------------------------------
;
;  Updated 15-Jul-2013
;	Make release independant
;


;--- Include MAKEMSI support (with my customisations and MSI branding) ------
#define VER_FILENAME.VER  Hercules-W32.Ver     
#define COMPANY_PRODUCT_ICON hercules.ico
#define COMPANY_SUMMARY_TEMPLATE Intel;1033
#define COMPANY_PROPERTY_MANUFACTURER  Hercules
#define COMPANY_WANT_TO_INSTALL_DOCUMENTATION   N
#include "HERCULES.MMH"


;--- Want to debug (not common) ---------------------------------------------
;#debug on
;#Option DebugLevel=^NONE, +OpSys^


;--- Define default location where file should install and add files --------
<$DirectoryTree Key="INSTALLDIR" Dir="[ProgramFilesFolder]\Hercules\R<$ProductVersion.MAJ>.<$ProductVersion.MIN>" CHANGE="\" PrimaryFolder="Y" Make="Y" Remove="y">

<$Files "..\msvc.DLLMOD.bin\*.dll" DestDir="INSTALLDIR">
<$Files "..\msvc.DLLMOD.bin\*.exe" DestDir="INSTALLDIR">
<$Files "..\html\*.*" SubDir='TREE' DestDir="[INSTALLDIR]\html\">
<$Files "..\util\*.*" DestDir="[INSTALLDIR]\util\">
<$Files "*.bat" DestDir="[INSTALLDIR]">
<$Files "..\README.*" DestDir="[INSTALLDIR]">
<$Files "..\release.notes" DestDir="[INSTALLDIR]">
<$Files "..\copyright" DestDir="[INSTALLDIR]">
<$Files "..\hercules.cnf" DestDir="[INSTALLDIR]">
<$Files "..\herclogo.txt" DestDir="[INSTALLDIR]">

<$msivalfilter ".D8D85FD0_537C_3A3A_9BEC_7A1B426637EC">

<$MergeModule "..\vc-redist\microsoft.vcxx.crt.x86_msm.msm" IgnoreErrors="Directory:SystemFolder">
 
;----------------------------------------------------------------------------
;--- Create Shortcut to "open" "setherc.bat" ---------------
;----------------------------------------------------------------------------
;<$Component "HercShortcuts" Create="Y" Directory_="INSTALLDIR" LM="Y" Make="y" Remove="y">

<$Component "HercShortcuts" Create="Y"  Directory_="INSTALLDIR" LM="Y">


   <$DirectoryTree Key="MENUDIR" Dir="[ProgramMenuFolder]\Hercules" Make="Y" Remove="y">
    
   <$DirectoryTree Key="XMENUDIRX" Dir="[ProgramMenuFolder]\Hercules\Hercules <$ProductVersion.MAJ>.<$ProductVersion.MIN> (32 Bit)" Make="Y" Remove="y">


    ;--- Add the shortcuts ----
    #(
        ;--- Shortcut in "Start->programs" -----------------------
        <$Shortcut
           Target="[INSTALLDIR]\herchere.bat"
           Arguments=^"[INSTALLDIR]"^
	   Dir="[XMENUDIRX]"  ;;Add to "Start->Programs
           Title="Hercules CMD Prompt (32 bit)"
           Description=^Open a Command Prompt with the Hercules Directory added to the path^
           Icon="@.\hercules.ico"
           WorkDir="INSTALLDIR"
	   Feature="."
        >
    #)


    ;--- Add the shortcut ----------------
    #(
        ;--- Shortcut in "Start->programs" -----------------------
        <$Shortcut
           Target="[INSTALLDIR]\html\index.html"
           Dir="[XMENUDIRX]"          ;;Add to "Start->Programs
           Title="Hercules Documentation"
           Description=^Open the Documentation in your default browser^
           WorkDir="INSTALLDIR"
           Icon="@.\rc_html.ico"
	   Feature="."
        >
    #)
<$/Component>