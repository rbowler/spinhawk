Building Hercules MSI files with Make MSI
=========================================

1. Download and install "makemsi" by following the instructions here:-

http://makemsi-manual.dennisbareis.com/install.htm

2. Update the ".ver" files as appropriate by ADDING the new version at the top.
   Use Start-Makemsi-GUID Generation to generate a new GUID and paste it
   into the "Hercules-W32.GUIDS" file. Then generate a second GUID and
   paste it into the "Hercules-W64.GUIDS" file.

3. Create a new folder in the root directory called "VC-redist" and copy in the Visual C redistributable files.
   vcredist_x86.exe
   vcredist_x64.exe

   The location of these depends on which Visual C you are using. For the "Microsoft Windows SDK V7.0" they are here

C:\Program Files\Microsoft SDKs\Windows\v7.0\Redist\VC

   You also need these merge modules:
   Microsoft_VC100_CRT_x86.msm
   Microsoft_VC100_CRT_x64.msm

   The merge modules should be copied into the VC-redist directory and renamed as:
   microsoft.vcxx.crt.x86_msm.msm
   microsoft.vcxx.crt.x64_msm.msm

   The merge modules for VS2010 can be downloaded from
   http://downloads.hercules-390.eu/vc2010-redist.zip

4. Download the compression and regular expression packages from:
   http://downloads.hercules-390.eu/zlib-1.2.5.zip
   http://downloads.hercules-390.eu/bzip2-1.0.6.zip
   http://downloads.hercules-390.eu/pcre-8.20.zip

   Unzip these packages to any suitable location (such as C:\Packages) and
   set environment variables in Computer-Properties-Advanced System Settings:
   SET ZLIB_DIR=C:\Packages\zlib-1.2.5
   SET BZIP2_DIR=C:\Packages\bzip2-1.0.6
   SET PCRE_DIR=C:\Packages\pcre-8.20

5. Open a Visual Studio Command Prompt and build the 32-bit version
   of Hercules in msvc.dllmod.bin:
   copy makefile.msvc makefile
   nmake clean
   nmake HERCVER=3.13

6. Open a Windows SDK Command Prompt and build the 64-bit version
   of Hercules in msvc.AMD64.bin:
   nmake clean
   nmake HERCVER=3.13

7. Browse to the Hercules-W32.MM or Hercules-W64.MM file
   and choose "Build MSI - Production"
