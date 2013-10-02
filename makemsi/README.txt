Building Hercules MSI files with Make MSI
=========================================

1. Download and install "makemsi" by following the instructions here:-

http://makemsi-manual.dennisbareis.com/install.htm

2. Update the ".ver" files as appropriate by ADDING the new version at the top.

3. Create a new folder in the root directory called "VC-redist" and copy in the Visual C redistributable files.
   The location of these depends on which Visual C you are using. For the "Microsoft Windows SDK V7.0" they are here

C:\Program Files\Microsoft SDKs\Windows\v7.0\Redist\VC

4. Browse to the Hercules-W32.mm or Hercules-W64.mm file and choose "Build MSI -Production"