#!/bin/bash

# Check if input file is provided and ends with .md
if [[ $# -ne 1 || "$1" != *.md ]]; then
  echo "Usage: $0 <filename>.md"
  exit 1
fi

# Extract the base filename (without extension)
input_file="$1"
base_name="${input_file%.md}"

# Output filename with .tex extension
output_file="${base_name}.tex"

# Run Pandoc to convert Markdown to LaTeX
pandoc -f markdown+tex_math_dollars -t latex -s -o "$output_file" "$input_file" --template=custom-template.tex --mathjax

# Check if Pandoc succeeded
if [[ $? -eq 0 ]]; then
  echo "Conversion successful: $output_file"
else
  echo "Pandoc conversion failed."
  exit 2
fi

