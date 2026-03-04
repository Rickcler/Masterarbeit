param(
    [string]$filePath = 'c:\Users\user\OneDrive\Desktop\Masterarbeit\Draft\Draft.tex'
)

# Read the file with proper encoding
$content = [System.IO.File]::ReadAllText($filePath, [System.Text.Encoding]::UTF8)

# Replace tabs with 4 spaces throughout the file
$content = $content -replace "`t", '    '

# Split into lines
$lines = $content -split "`n"
$formattedLines = @()

foreach ($line in $lines) {
    # Remove trailing whitespace
    $line = $line -replace '\s+$', ''
    $formattedLines += $line
}

# Join back with consistent newlines
$formatted = $formattedLines -join "`n"

# Ensure file ends with single newline
$formatted = $formatted.TrimEnd() + "`n"

# Save with UTF-8 encoding
[System.IO.File]::WriteAllText($filePath, $formatted, [System.Text.Encoding]::UTF8)
Write-Host "Formatting complete:"
Write-Host "- Tabs replaced with 4 spaces"
Write-Host "- Trailing whitespace removed"
Write-Host "- Consistent line endings"
