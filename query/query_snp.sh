if [[ $# -eq 0 ]] ; then
    printf "Querying the ENTEx AS catalog.\nPlease enter the path to the input file, followed by at most three of the following types of keywords:\n      Assay (e.g. H3K27ac)\n      Individual (e.g. ENC-001)\n      Tissue (e.g. spleen)\nThe output is all AS SNPs in the catalog from the query.\n"
    exit 0
fi

f="$1"
q1="$2"
q2="$3"
q3="$4"


echo "Querying $q1|$q2|$q3 on file $f"
grep "$q1" $f --color=never | grep "$q2" --color=never | grep "$q3" --color=never | awk '{print $1,$2,$4,$NF}' | sort | uniq |  awk '{if($NF==1) print $0}'
