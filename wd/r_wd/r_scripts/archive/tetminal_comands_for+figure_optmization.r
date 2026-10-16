gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dCompatibilityLevel=1.5 -dPDFSETTINGS=/screen -sOutputFile=output.pdf input.pdf



gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dCompatibilityLevel=1.5 -dPDFSETTINGS=/screen -sOutputFile=for_CLAUDE_test_horizontal_eye_legend.pdf 

# Quality/Compression Levels ('-dPDFSETTINGS')
# /screen   - lowest quality, smallest file size (72 dpi)
# /ebook    - medium quality (150 dpi)
# /printer  - higher quality (300 dpi)
# /prepress - highest quality, large file size (300 dpi)

# Additional Compression Options
-dColorConversionStrategy=/RGB   # Convert colors to RGB
-dColorImageDownsampleType=/Bicubic  # Best downsampling method
-dColorImageResolution=150       # Set color image resolution
-dGrayImageDownsampleType=/Bicubic
-dGrayImageResolution=150
-dMonoImageDownsampleType=/Bicubic
-dMonoImageResolution=300

# Embed Fonts
-dEmbedAllFonts=true
-dSubsetFonts=true
-dMaxSubsetPct=100