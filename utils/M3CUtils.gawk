##
# @brief Convert the M3CLabel like this "C_1+C_1+C_2"
#        to its reduced version "2C_1/C_2"
#
function reduceM3CLabel( M3CLabel )
{
	split(M3CLabel,arr,"+")
	
	for( i=1; i<=length(arr); i++ )
		if( arr[i] in map )
			map[arr[i]]+=1
		else
			map[arr[i]]=1
	
	i=1
	for( item in map ){
		keys[i]=item
		i++
	}
	
	n = asort(keys)
	
	newLabel=""
	for( i=1; i<=n; i++ ){
		if( map[keys[i]]!=1 )
			newLabel = newLabel""map[keys[i]]""keys[i]
		else
			newLabel = newLabel""keys[i]
			
		if( i!=n )
			newLabel = newLabel"/"
	};
	
	delete map
	delete arr
	delete keys
	
	return newLabel
}
