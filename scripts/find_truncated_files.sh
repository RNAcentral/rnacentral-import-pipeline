# find all ncr files that might be truncated
# requires full path as a parameter
for file in "$1"/*.ncr
do
	if tail -2 "$file" | grep -q ^\/\/
	then
	  :
	else
	  echo "No match with the pattern in $file"
	fi
done
