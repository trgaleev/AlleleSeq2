if [[ $# -eq 0 ]] ; then
    printf "Querying the hetSNV catalog.\nPlease enter the path to the input file, followed by at most three of the following types of keywords:\n      Assay (e.g. H3K27ac)\n      Individual (e.g. ENC-001)\n      Tissue (e.g. spleen)\nThe output will count the number of AS SNPs in the catalog from the query.\n"
    exit 0
fi

f="$1"
q1="$2"
q2="$3"
q3="$4"


echo "Querying $q1|$q2|$q3 on file $f"
grep "$q1" $f | grep "$q2" | grep "$q3" | awk '{print $1,$2,$4,$NF}' | sort | uniq | awk '{print $NF}' | sort | uniq  -c | awk '{if($NF==1) print "AS "$1}'
