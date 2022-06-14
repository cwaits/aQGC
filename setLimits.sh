#!/bin/bash -x
# run in aQGC dir
# sample naming convention to automatically pick correct normalization in analyzeDelphes.py: INT6_T1.root with 6 the CoM energy and T1 the aQGC operator
if [[ "$PWD" != "/raid01/users/cwaits/aQGC" ]];
then
    exit
fi
# Example: ./setLimits.sh -i $PWD/path/to/INT_sample -q $PWD/path/to/QUAD_sample -s $PWD/path/to/SM_sample  -t $PWD/example/dir/  -h R_dijet_mass
# -g YES  -c cut_dictionary"
# the tag is the directory that will hold the INT, QUAD, and background samples to be sent eft-fun
# YES for gridscanner mode will not ask for user input, only use when running this script
# through the separate gridscan scrpt
while getopts i:q:s:t:g:h:c: flag; do
    echo "flag -$flag, Argument $OPTARG";
    case "$flag" in
        i) INT_path=$OPTARG;;
        q) QUAD_path=$OPTARG;;
        s) SM_path=$OPTARG;;
        t) tag=$OPTARG;;
        g) gridscannerMode=$OPTARG;;
        h) target_hist=$OPTARG;;
        c) cut_dict=$OPTARG;;
    esac
done

INT_name="$(basename -- $INT_path)"
QUAD_name="$(basename -- $QUAD_path)"
SM_name="$(basename -- $SM_path)"

INT_norm="${INT_name/.root/}"
QUAD_norm="${QUAD_name/.root/}"
SM_norm="${SM_name/.root/}"

if [[ "$gridscannerMode" == "YES" ]];
#removes contents of $tag in gridscanner mode
then
    if [[ -d "$tag" ]];
    then
        rm -rf "$tag"
        mkdir "$tag"
    else
        mkdir "$tag"
    fi
    #run with cut dictionary when in gridscanner mode
    #python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $INT_path $INT_norm $cut_dict
    #python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $QUAD_path $QUAD_norm $cut_dict
    #python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $SM_path $SM_norm $cut_dict
else
    #remove contents of output directory before moving histogram files
    if [[ -d "$tag" ]];
    then
        echo "$tag"
        echo "Remove contents of above directory?"
        read -p "Enter YES or NO: " uservar
    fi

    if [[ "$uservar" == "YES" ]];
    then
        rm -rf "$tag"
        mkdir "$tag"
    else
        exit
    fi
    #make output directory if it does not exist
    if [[ ! -d "$tag" ]];
    then
        mkdir "$tag"
    fi
    python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $INT_path $INT_norm
    python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $QUAD_path $QUAD_norm
    python2.7 /raid01/users/cwaits/aQGC/analyzeDelphes.py $SM_path $SM_norm
fi

# cp analyzeDelphes.py output to tag dir
INT_path="${INT_path/.root/.hist.root}"
QUAD_path="${QUAD_path/.root/.hist.root}"
SM_path="${SM_path/.root/.hist.root}"

INT_name="$(basename -- $INT_path)"
QUAD_name="$(basename -- $QUAD_path)"
SM_name="$(basename -- $SM_path)"

cp "$INT_path" "$tag"
cp "$QUAD_path" "$tag"
cp "$SM_path" "$tag"

INT_path= "$tag"/"$INT_name"
QUAD_path="$tag"/"$QUAD_name"
SM_path="$tag"/"$SM_name"

# merge low content bins
python2.7 binMerger.py "$SM_path" "$target_hist"

python2.7 reBin.py "$SM_path" "$INT_path" "$target_hist"

python2.7 reBin.py "$SM_path" "$QUAD_path" "$target_hist"

# edit eft-fun config and set limits
cd ../eft-fun
var1="histopath = "
var2="$var1$tag"
echo "$var2"
sed -i "3 d" configs/mumu_vbs/mumu_vbs_mumu_mjj.cfg
sed -i "3 i $var2" configs/mumu_vbs/mumu_vbs_mumu_mjj.cfg
./bin/eftfun.py -m scan -a -p T1 -c all -i configs/mumu_vbs/mumu_vbs_mumu.cfg
cd ../aQGC 
