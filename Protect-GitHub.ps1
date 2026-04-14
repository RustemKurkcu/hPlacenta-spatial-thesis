Write-Host "Scanning Placenta Pipeline for massive files (>100MB)..." -ForegroundColor Cyan

# Find files larger than 100MB, ignoring the hidden Git folder
$largeFiles = Get-ChildItem -File -Recurse | Where-Object { $_.Length -ge 100MB -and $_.FullName -notmatch "\\.git\\" }

if ($largeFiles) {
    Write-Host "Found $($largeFiles.Count) files over the 100MB limit. Blocking them now..." -ForegroundColor Yellow
    foreach ($file in $largeFiles) {
        # Convert the file path to a Git-friendly format
        $relativePath = $file.FullName.Replace($PWD.Path + "\", "").Replace("\", "/")
        
        # Check if the file is already blocked in .gitignore
        $ignoreContent = If (Test-Path .gitignore) { Get-Content .gitignore } Else { @() }
        if ($ignoreContent -notcontains $relativePath) {
            Add-Content -Path .gitignore -Value $relativePath
            Write-Host "BLOCKED: $relativePath" -ForegroundColor Red
        }
    }
} else {
    Write-Host "All files are under 100MB! You are safe to push." -ForegroundColor Green
}

Write-Host "Done! You can now run your git add, commit, and push commands." -ForegroundColor Cyan
