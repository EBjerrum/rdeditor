for file in *.png; do
    convert "$file" -negate "$file"
done