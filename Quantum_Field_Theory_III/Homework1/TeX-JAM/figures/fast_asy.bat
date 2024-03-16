@echo off

set "folder_path=C:\Users\pgraham1\Documents\GitHub\PSI\Classical_Physics\Homework1\TeX-JAM\figures"

for %%f in ("%folder_path%\*.asy") do (
    echo Running command on file: %%f
    asy -noV %%f
)

for %%f in ("%folder_path%\*.eps") do (
    echo converting to pdf: %%f
    epstopdf %%f
)





