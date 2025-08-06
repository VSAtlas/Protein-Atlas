@echo off
REM ---- Set up Conda environment for protein_automation ----

REM Check if environment exists
CALL conda info --envs | findstr "docking-env" >nul
IF ERRORLEVEL 1 (
    echo Creating new environment: docking-env
    CALL conda create -n docking-env python=3.10 -y
) ELSE (
    echo Environment 'docking-env' already exists.
)

REM Activate the environment
CALL conda activate docking-env

REM Confirm interpreter path
echo Active Python Interpreter:
where python

REM Install conda packages from conda-forge
echo Installing conda packages (rdkit, openbabel, tqdm)...
CALL conda install -c conda-forge rdkit openbabel tqdm -y
IF ERRORLEVEL 1 (
    echo Conda package installation failed. Exiting.
    pause
    exit /b 1
)

REM Install pip packages
echo Installing pip packages (pandas, biopython, p2ranktools, pywin32)...
CALL pip install pandas biopython p2ranktools pywin32
IF ERRORLEVEL 1 (
    echo Pip package installation failed. Exiting.
    pause
    exit /b 1
)

REM Verify Biopython install
echo Verifying Biopython installation...
python -c "from Bio import SeqIO; print('Biopython is installed and working.')" || (
    echo Biopython test failed. Exiting.
    pause
    exit /b 1
)

echo.
echo Environment setup complete.
echo To activate later: conda activate docking-env
pause
