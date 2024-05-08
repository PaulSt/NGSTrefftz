$ErrorActionPreference = "Stop"

if (test-path _skbuild) {
    cmd.exe /c rd /s /q _skbuild
}
if (test-path dist) {
    cmd.exe /c rd /s /q dist
}
if (test-path venv_ngs) {
    cmd.exe /c rd /s /q venv_ngs
}

##python --version

#$py=$args[0]
#& $py\python.exe -m venv venv_ngs
#venv_ngs\scripts\Activate.ps1
#$env:PATH += ";$env:CI_PROJECT_DIR\venv\bin"
#$env:NGSolve_DIR = Join-Path (Get-Item .).FullName "venv_ngs\Lib\site-packages\ngsolve\cmake"
#$env:Netgen_DIR = Join-Path (Get-Item .).FullName "venv_ngs\Lib\site-packages\netgen\cmake"
#$env:LIB = "$env:LIB;" + (Join-Path (Get-Item .).FullName "venv_ngs\Lib")
#$env:LIBPATH = "$env:LIBPATH;" + (Join-Path (Get-Item .).FullName "venv_ngs\Lib")
#$env:PYDIR = (Join-Path (Get-Item .).FullName "venv_ngs\")
#Get-ChildItem venv_ngs
#Get-ChildItem venv_ngs\lib
#Get-ChildItem venv_ngs\include
#Get-ChildItem venv_ngs\scripts
#Write-Host $env:LIB


pip3 install scikit-build wheel numpy twine mkl-devel==2022.* mkl==2022.* setuptools setuptools_scm
pip3 install -r ngsolve_version.txt


#$env:PYDIR = "$env:Python3_ROOT_DIR"
#$env:PATH += ";$env:pythonLibrary;$env:Python3_ROOT_DIR"
$env:NGSolve_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\ngsolve\cmake"
$env:Netgen_DIR = "$env:Python3_ROOT_DIR\lib\site-packages\netgen\cmake"
#$env:LIB = "$env:LIB;$env:Python3_ROOT_DIR\libs"
#$env:LIBPATH = "$env:LIBPATH;$env:Python3_ROOT_DIR\libs"

#Get-ChildItem $env:pythonLibrary
#Get-ChildItem $env:Python3_ROOT_DIR
#Get-ChildItem $env:Python3_ROOT_DIR\lib
#Get-ChildItem $env:Python3_ROOT_DIR\libs
#Get-ChildItem $env:Python3_ROOT_DIR\include
#Get-ChildItem $env:Python3_ROOT_DIR\scripts
#Get-ChildItem $env:Python3_ROOT_DIR\share
#Get-ChildItem $env:Python3_ROOT_DIR\lib\site-packages
#Get-ChildItem $env:Python3_ROOT_DIR\lib\site-packages\ngsolve

Set-Location ../..
Get-ChildItem
python setup.py bdist_wheel -G"Visual Studio 16 2019" -d wheelhouse
