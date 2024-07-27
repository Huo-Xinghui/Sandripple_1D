set OPTC=-I"C:\Program Files (x86)\Microsoft SDKs\MPI\Include\intel"
set OPTL="C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\msmpi.lib" "C:\Program Files (x86)\Microsoft SDKs\MPI\Lib\x64\msmpifec.lib"
call ifx %OPTC% -o %1.exe %1.f90 %OPTL%