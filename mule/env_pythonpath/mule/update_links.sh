#! /bin/bash


function create_links () {
	# create_links [level_id] [file1] [file2] [file3]
	local LEVEL=$1
	for i in ${@:2}; do
		if [[ "$i" == *__pycache__ ]]; then
			continue
		fi
		if [[ "$i" == *mule_local/python/__init__.py ]]; then
			continue
		fi

		if [[ "$LEVEL" != "0" ]]; then
			if [[ -e "$i" ]]; then
				echo "FILE $i"
				ln "${i}" ./ -sf || exit 1
				continue
			fi
		else
			local NAME=${i##*/}
			if [[ -d "$i" ]]; then
				echo "DIR $i"

				mkdir -p "$NAME" || exit 1

				cd "$NAME" || exit 1
				create_links $((LEVEL+1)) ../$i/* || exit 1

				cd "../" || exit 1
				continue
			fi

			if [[ "$i" == *.py ]]; then
				echo "FILE $i"
				ln "${i}" ./ -sf || exit 1
				continue
			fi
		fi

		echo "Don't know how to process $i"
		exit 1
	done
}

for i in *; do
	# Remove all directories
	if [[ -d "$i" ]]; then
		rm -r "$i"
	fi
	# Remove all symlinks
	if [[ -h "$i" ]]; then
		rm -r "$i"
	fi
done

create_links 0 ../../python/*

create_links 0 ../../../mule_local/python/*

