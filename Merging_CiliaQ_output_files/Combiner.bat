@echo off
setlocal

set "firstfile=true"
for %%F in (*CQs.txt) do (
    if defined firstfile (
        type "%%F" > AllCQsFilesMerged.txt
        set "firstfile="
    ) else (
        echo. >> AllCQsFilesMerged.txt
        type "%%F" >> AllCQsFilesMerged.txt
    )
)