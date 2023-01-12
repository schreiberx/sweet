#! /bin/bash


function create_links () {
	for i in $@; do
		if [[ "$i" == *__pycache__ ]]; then
			continue
		fi
		if [[ "$i" == *__init__.py ]]; then
			continue
		fi

		if [[ -e "$i" ]]; then
			echo "$i"
			ln "${i}" ./ -sf
		fi
	done
}


create_links ../../python/*

create_links ../../../mule_local/python/*

