param(
    [string]$filePath = 'c:\Users\user\OneDrive\Desktop\Masterarbeit\Draft\Draft.tex'
)

# Read the file
$content = Get-Content -Path $filePath -Raw

# Replace tabs with 4 spaces
$content = $content -replace "`t", '    '

# Normalize multiple consecutive spaces at line start (but keep logical indentation)
# This regex finds lines and normalizes their leading whitespace
$lines = $content -split "`n"
$formattedLines = @()

foreach ($line in $lines) {
    # Remove trailing whitespace
    $line = $line -replace '\s+$', ''
    
    # Count leading spaces and convert to multiples of 4
    if ($line -match '^(\s+)') {
        $spaces = $matches[1].Length
        # Round to nearest multiple of 4
        $indentLevel = [Math]::Round($spaces / 4)
        $indentation = '    ' * $indentLevel
        $trimmed = $line -replace '^\s+', ''
        $formattedLines += $indentation + $trimmed
    } else {
        $formattedLines += $line
    }
}

# Join lines back and ensure proper newlines
$formatted = $formattedLines -join "`n"

# Ensure file ends with a single newline
$formatted = $formatted -replace '\n+$', "`n"

# Save back
[System.IO.File]::WriteAllText($filePath, $formatted, [System.Text.Encoding]::UTF8)
Write-Host "File formatted successfully: $filePath"
Write-Host "- Tabs replaced with spaces"
Write-Host "- Indentation normalized to 4-space increments"
Write-Host "- Trailing whitespace removed"
