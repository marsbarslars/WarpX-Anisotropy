@echo off
REM Windows wrapper around warpx.py. The driver is stdlib-only, so the system
REM Python works — .venv does not exist yet when `warpx sync` first runs.
REM
REM Note the goto structure: %ERRORLEVEL% inside a parenthesised block expands
REM when the block is parsed, not when it runs, so it would report a stale exit
REM code. Keeping `exit /b` out of any block is what makes it accurate.
setlocal
where py >nul 2>&1 && goto usepy
where python >nul 2>&1 && goto usepython
echo warpx: no Python 3.8+ interpreter found on PATH 1>&2
exit /b 1

:usepy
py -3 "%~dp0warpx.py" %*
exit /b %ERRORLEVEL%

:usepython
python "%~dp0warpx.py" %*
exit /b %ERRORLEVEL%
