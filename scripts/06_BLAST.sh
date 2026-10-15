#!/bin/bash
#SBATCH --job-name=BLASTn
#SBATCH --partition=medium
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32

# =============================================================================
# BLASTn
#
# Each query sequence is BLASTed against the DB, and the input is split into
# <basename>_hits.fasta and <basename>_no_hits.fasta. Alignment gaps are removed for
# the search, but the original sequences (with gaps) are written to the output.
#
# Queries are batched into chunks, and the chunks run concurrently. The
# database is therefore scanned once per chunk instead of once per sequence.
# The number of chunks adapts to the input so that all allocated CPUs are used.
#
# - BLAST sees synthetic query IDs (q1, q2, ...), which are mapped back to the
#   real IDs afterwards. Duplicate or truncated headers cannot collide.
# - Per-sequence .out files are written atomically (.part, then renamed). An
#   empty .out file therefore always means "no hits".
# - If a chunk fails, its sequences are retried one at a time. Sequences that
#   still fail are listed in <basename>_failed.tsv.
# - Sequences that already have a .out file are not BLASTed again,
#   so a re-run only reruns failures.
#
# Usage:
#   sbatch blastn_screen.sh -i INPUT -o OUTPUT_DIR -d DB [-n DB_NAME]
#   bash   blastn_screen.sh -i INPUT -o OUTPUT_DIR -d DB [-n DB_NAME]
#
#   -i, --input     Directory of FASTA files (.fa/.fasta/.fas) or a single multi-FASTA file
#   -o, --output    Output directory
#   -d, --db        Directory containing a prebuilt BLAST nucleotide database,
#                   or a multi-FASTA file to build a database from
#   -n, --db-name           Database name, without file extensions (optional, see below)
#   -t, --taxids            Optional. Only search these taxa (comma-separated NCBI taxids)
#   -x, --negative-taxids   Optional. Search everything EXCEPT these taxa (comma-separated)
#   -h, --help              Show usage
#
# The same settings can be set in the CONFIGURATION block below instead.
# A flag on the command line overrides the value in the block.
#
# Taxonomy filtering (-t / -x):
#   - Only one of the two can be used. I.e. they are mutually exclusive.
#   - Needs a prebuilt NCBI (v5) database with taxonomy built in (e.g. core_nt)
#     Not available when -d is a FASTA file.
#   - taxonomy4blast.sqlite3 must be in the database directory or in $BLASTDB.
#     It ships with the NCBI databases and in taxdb.tar.gz.
#   - If BLAST rejects the taxids (e.g. a taxid that is not in the database),
#     the script stops instead of retrying every sequence individually.
#
# If DB_NAME is not set:
#   - DB is a FASTA file: the FASTA filename without extension is used.
#   - DB is a directory:  the database is detected automatically if the
#     directory contains exactly one; otherwise the script stops and lists them.
#
# Outputs per input file <basename>:
#   OUTPUT_DIR/<basename>/<sanitised_header>.out  per-sequence BLAST tabular (outfmt 6)
#   OUTPUT_DIR/<basename>_blast.tsv               all hits combined, with a header row
#   OUTPUT_DIR/<basename>_hits.fasta              sequences with >=1 hit
#   OUTPUT_DIR/<basename>_no_hits.fasta           sequences with no hits
#   OUTPUT_DIR/<basename>_failed.tsv              only written if something failed
#
# =============================================================================

set -o pipefail

# Record start time
start_time=$(date +%s)

# ======= CONFIGURATION =======
# The next four can also be given on the command line (-i / -o / -d / -n),
# which overrides the values here. Leave empty if they are always passed as flags.

# Input: a directory of FASTA files (.fa/.fasta/.fas) or a single multi-FASTA file
INPUT_PATH=""

# Output directory where BLAST results will be saved
OUTPUT_DIR=""

# Path to the NCBI database directory, or a multi-FASTA file to build a database from
DB_DIR=""

# Name of the BLAST database to use (without file extensions).
# Leave empty to determine it automatically (see header).
# If DB_DIR is a FASTA file, this will be the name of the created database.
DB_NAME=""

# Optional taxonomy filter (comma-separated NCBI taxids, e.g. "561,590").
# Set at most one of the two. Overridden by -t / -x. Leave both empty for no filter.
TAXIDS=""            # only search these taxa
NEGATIVE_TAXIDS=""   # search everything except these taxa

# Type of BLAST database to create (nucl or prot); only used if DB_DIR is a FASTA file
DB_TYPE="nucl"

# Additional BLAST options. Do NOT put -num_threads here: the script sets the
# thread count for each blastn process itself.
BLAST_OPTIONS="-evalue 1e-5 -max_target_seqs 50"

# Total CPUs to use (defaults to the SLURM allocation)
THREADS="${SLURM_CPUS_PER_TASK:-32}"

# Upper limit on query sequences per blastn call. Normally the input is split
# into THREADS chunks; this cap only matters for very large inputs.
MAX_CHUNK_SIZE=500
# ===================================================

# -----------------------------------------------------------------------------
# Command-line arguments
# -----------------------------------------------------------------------------
usage() {
    cat <<EOF
Usage: $(basename "$0") -i INPUT -o OUTPUT_DIR -d DB [-n DB_NAME]

  -i, --input     Directory of FASTA files (.fa/.fasta/.fas) or a single multi-FASTA file
  -o, --output    Output directory
  -d, --db        Directory containing a prebuilt BLAST nucleotide database,
                  or a multi-FASTA file to build a database from
  -n, --db-name   Database name, without file extensions. Optional:
                    DB is a FASTA     -> defaults to the FASTA filename without extension
                    DB is a directory -> auto-detected if it contains exactly one database
  -t, --taxids           Optional. Only search these taxa (comma-separated NCBI taxids)
  -x, --negative-taxids  Optional. Search everything except these taxa (comma-separated)
                         -t and -x cannot be used together. Both need a prebuilt NCBI
                         database with taxonomy (e.g. core_nt), not a FASTA -d.
  -h, --help      Show this message

Flags override the values in the script's CONFIGURATION block.
A flag given an empty value (e.g. -n "") is treated as not given.
EOF
}

# need_value FLAG N_REMAINING_ARGS VALUE
need_value() {
    if (( $2 < 2 )) || [[ "$3" == -* ]]; then
        echo "Error: $1 requires a value" >&2
        usage >&2
        exit 1
    fi
}

while (( $# > 0 )); do
    case "$1" in
        -i|--input)   need_value "$1" $# "$2"; [ -n "$2" ] && INPUT_PATH="$2"; shift 2 ;;
        -o|--output)  need_value "$1" $# "$2"; [ -n "$2" ] && OUTPUT_DIR="$2"; shift 2 ;;
        -d|--db)      need_value "$1" $# "$2"; [ -n "$2" ] && DB_DIR="$2";     shift 2 ;;
        -n|--db-name) need_value "$1" $# "$2"; [ -n "$2" ] && DB_NAME="$2";    shift 2 ;;
        -t|--taxids)  need_value "$1" $# "$2"; [ -n "$2" ] && TAXIDS="$2";     shift 2 ;;
        -x|--negative-taxids)
                      need_value "$1" $# "$2"; [ -n "$2" ] && NEGATIVE_TAXIDS="$2"; shift 2 ;;
        -h|--help)    usage; exit 0 ;;
        *)
            echo "Error: unknown argument '$1'" >&2
            usage >&2
            exit 1
            ;;
    esac
done

missing=()
[ -z "$INPUT_PATH" ] && missing+=("-i/--input  (INPUT_PATH)")
[ -z "$OUTPUT_DIR" ] && missing+=("-o/--output (OUTPUT_DIR)")
[ -z "$DB_DIR" ]     && missing+=("-d/--db     (DB_DIR)")
if (( ${#missing[@]} > 0 )); then
    echo "Error: missing required setting(s), set them as flags or in the CONFIGURATION block:" >&2
    printf '  %s\n' "${missing[@]}" >&2
    echo >&2
    usage >&2
    exit 1
fi

# Taxonomy filter: at most one of TAXIDS / NEGATIVE_TAXIDS, as comma-separated integers
TAXIDS="${TAXIDS//[[:space:]]/}"
NEGATIVE_TAXIDS="${NEGATIVE_TAXIDS//[[:space:]]/}"
TAX_ARGS=()
TAX_DESC="none"
if [ -n "$TAXIDS" ] && [ -n "$NEGATIVE_TAXIDS" ]; then
    echo "Error: --taxids and --negative-taxids cannot be used together (BLAST does not allow it)." >&2
    echo "  TAXIDS=$TAXIDS" >&2
    echo "  NEGATIVE_TAXIDS=$NEGATIVE_TAXIDS" >&2
    echo "Check both the command line and the CONFIGURATION block." >&2
    exit 1
fi
taxid_re='^[0-9]+(,[0-9]+)*$'
if [ -n "$TAXIDS" ]; then
    if [[ ! "$TAXIDS" =~ $taxid_re ]]; then
        echo "Error: --taxids must be comma-separated numeric taxids (e.g. 561,590), got: $TAXIDS" >&2
        exit 1
    fi
    TAX_ARGS=(-taxids "$TAXIDS")
    TAX_DESC="only taxids $TAXIDS"
elif [ -n "$NEGATIVE_TAXIDS" ]; then
    if [[ ! "$NEGATIVE_TAXIDS" =~ $taxid_re ]]; then
        echo "Error: --negative-taxids must be comma-separated numeric taxids (e.g. 561,590), got: $NEGATIVE_TAXIDS" >&2
        exit 1
    fi
    TAX_ARGS=(-negative_taxids "$NEGATIVE_TAXIDS")
    TAX_DESC="excluding taxids $NEGATIVE_TAXIDS"
fi

# A taxonomy filter in BLAST_OPTIONS as well would conflict with the flags
tax_opt_re='(^|[[:space:]])-(negative_)?taxid(s|list)([[:space:]]|$)'
if (( ${#TAX_ARGS[@]} > 0 )) && [[ "$BLAST_OPTIONS" =~ $tax_opt_re ]]; then
    echo "Error: BLAST_OPTIONS already contains a taxonomy filter (-taxids/-taxidlist/-negative_taxids/-negative_taxidlist)." >&2
    echo "Remove it from BLAST_OPTIONS or don't use --taxids/--negative-taxids." >&2
    exit 1
fi

# Set path explicitly
export PATH="/mnt/apps/users/dparsons/conda/bin:$PATH"

# Source conda.sh & activate env
source /mnt/apps/users/dparsons/conda/etc/profile.d/conda.sh
conda activate blast

BLAST_COLS="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if blastn is available
if ! command -v blastn &> /dev/null; then
    echo "Error: blastn command not found. Please make sure BLAST+ is installed and in your PATH"
    exit 1
fi

# Remove -num_threads from BLAST_OPTIONS if present (threads are set per chunk)
if [[ "$BLAST_OPTIONS" =~ -num_threads ]]; then
    echo "Warning: removing -num_threads from BLAST_OPTIONS (threads are set per chunk by this script)"
    BLAST_OPTIONS=$(sed -E 's/-num_threads[[:space:]]+[0-9]+//g' <<< "$BLAST_OPTIONS")
fi

# Check if database path exists
if [ ! -e "$DB_DIR" ]; then
    echo "Error: Database path $DB_DIR does not exist"
    exit 1
fi

# -----------------------------------------------------------------------------
# list_db_names DIR
#   Prints the name of each nucleotide BLAST database directly inside DIR.
#   Multi-volume databases (name.nal + name.00.nhr, name.01.nhr, ...) are
#   listed once.
# -----------------------------------------------------------------------------
list_db_names() {
    find "$1" -maxdepth 1 -type f \( -name '*.nal' -o -name '*.nhr' \) -printf '%f\n' \
        | sed -E 's/\.nal$//; s/\.[0-9]+\.nhr$//; s/\.nhr$//' \
        | sort -u
}

# Variables for database handling
DB_PATH=""

# Determine if DB_DIR is a file (FASTA) or a directory (prebuilt database)
if [ -f "$DB_DIR" ]; then
    # Check if the file looks like a FASTA file (starts with ">")
    if grep -q "^>" "$DB_DIR"; then
        echo "Detected that DB_DIR is a FASTA file."

        # Databases built here have no taxonomy, so a taxid filter cannot work
        if (( ${#TAX_ARGS[@]} > 0 )); then
            echo "Error: --taxids/--negative-taxids need a prebuilt NCBI database with taxonomy (e.g. core_nt)."
            echo "       DB_DIR is a FASTA file, and the database built from it has no taxonomy information."
            exit 1
        fi

        # Check if makeblastdb is available
        if ! command -v makeblastdb &> /dev/null; then
            echo "Error: makeblastdb command not found but needed to create a database. Please make sure BLAST+ is installed and in your PATH"
            exit 1
        fi

        # Get the directory and base name of the FASTA file
        fasta_dir=$(dirname "$DB_DIR")
        fasta_basename=$(basename "$DB_DIR")
        fasta_name="${fasta_basename%.*}"

        # If DB_NAME is not set, use the FASTA filename instead
        if [ -z "$DB_NAME" ]; then
            DB_NAME="$fasta_name"
            echo "DB_NAME not set; using FASTA filename: $DB_NAME"
        fi

        # Put the database in the same directory as the FASTA file
        DB_PATH="${fasta_dir}/${DB_NAME}"

        # Check if the database already exists
        if [ -e "${DB_PATH}.nhr" ] || [ -e "${DB_PATH}.00.nhr" ] || [ -e "${DB_PATH}.nal" ]; then
            echo "BLAST database already exists at $DB_PATH. Using existing database."
        else
            echo "Creating BLAST database from $DB_DIR in the same directory..."
            if ! makeblastdb -in "$DB_DIR" -dbtype "$DB_TYPE" -out "$DB_PATH"; then
                echo "Error: Failed to create BLAST database from $DB_DIR"
                exit 1
            fi
            echo "BLAST database created successfully at $DB_PATH"
        fi
    else
        echo "Error: The file at $DB_DIR does not appear to be a FASTA file (missing '>' header)"
        exit 1
    fi
else
    # It's a directory, so assume it contains a prebuilt database

    # If DB_NAME is not set, use the database in the directory if there is exactly one
    if [ -z "$DB_NAME" ]; then
        mapfile -t db_names < <(list_db_names "$DB_DIR")
        if (( ${#db_names[@]} == 1 )); then
            DB_NAME="${db_names[0]}"
            echo "DB_NAME not set; auto-detected database: $DB_NAME"
        elif (( ${#db_names[@]} == 0 )); then
            echo "Error: DB_NAME not set and no BLAST nucleotide database found in $DB_DIR"
            exit 1
        else
            echo "Error: DB_NAME not set and $DB_DIR contains ${#db_names[@]} databases."
            echo "Set DB_NAME or pass -n with one of:"
            printf '  %s\n' "${db_names[@]}"
            exit 1
        fi
    fi

    DB_PATH="${DB_DIR}/${DB_NAME}"

    # Check the database exists by looking for common database file extensions
    if [ ! -e "${DB_PATH}.nhr" ] && [ ! -e "${DB_PATH}.00.nhr" ] && [ ! -e "${DB_PATH}.nal" ]; then
        echo "Error: BLAST database '${DB_NAME}' not found at $DB_DIR"
        echo "Available databases in $DB_DIR:"
        list_db_names "$DB_DIR" | sed 's/^/  /'
        exit 1
    fi

    # Without taxonomy4blast.sqlite3, BLAST can warn and then search the whole
    # database unfiltered, so require it up front.
    if (( ${#TAX_ARGS[@]} > 0 )); then
        tax_sqlite=""
        tax_dirs=( "$DB_DIR" )
        if [ -n "${BLASTDB:-}" ]; then
            IFS=: read -r -a blastdb_dirs <<< "$BLASTDB"
            tax_dirs+=( "${blastdb_dirs[@]}" )
        fi
        for d in "${tax_dirs[@]}"; do
            if [ -n "$d" ] && [ -f "$d/taxonomy4blast.sqlite3" ]; then
                tax_sqlite="$d/taxonomy4blast.sqlite3"
                break
            fi
        done
        if [ -z "$tax_sqlite" ]; then
            echo "Error: taxonomy filter requested but taxonomy4blast.sqlite3 was not found in:"
            printf '  %s\n' "${tax_dirs[@]}"
            echo "It ships with the NCBI preformatted databases and in"
            echo "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz - put it alongside the database."
            exit 1
        fi
        echo "Taxonomy filter: $TAX_DESC (using $tax_sqlite)"
    fi
fi

# Scratch space for chunk files, removed on exit. Per-sequence results are
# written straight to OUTPUT_DIR, so nothing needed for resuming lives here.
WORK_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/blast_screen.XXXXXX")
cleanup() {
    local pids
    pids=$(jobs -p)
    [ -n "$pids" ] && kill $pids 2>/dev/null
    rm -rf "$WORK_ROOT"
}
trap cleanup EXIT
trap 'echo "Interrupted - completed sequences are kept; re-run to resume"; exit 130' INT TERM

# Created (with BLAST's error message) if BLAST rejects the taxonomy filter.
# Chunks check it so that no further work is started.
ABORT_FILE="$WORK_ROOT/abort_taxonomy"

# is_taxonomy_error ERRFILE
#   True if a taxonomy filter is in use and BLAST's error mentions taxonomy/taxids
#   (e.g. "Taxonomy ID(s) not found", or a v4 database without taxonomy).
is_taxonomy_error() {
    (( ${#TAX_ARGS[@]} > 0 )) && grep -qiE 'taxonom|taxid' "$1" 2>/dev/null
}


# -----------------------------------------------------------------------------
# scan_fasta
#   Reads the input once. For every record it writes a line to /map.tsv:
#     synthetic_id  qseqid  out_path  status  original_header
#   status: done    = .out already exists (resume)
#           empty   = no sequence left after removing gaps
#           pending = needs BLAST
#   Pending records are written to /pending.fa with gaps removed and
#   synthetic IDs as headers.
# -----------------------------------------------------------------------------
scan_fasta() {
    awk -v outdir="$2" -v map="$3/map.tsv" -v pend="$3/pending.fa" '
        function sanitize(h,   s) {
            s = h
            sub(/^>/, "", s)
            gsub(/[\/|:*"<>?[:space:]]/, "_", s)
            return substr(s, 1, 100)
        }
        function flush(   r, tmp, status) {
            if (n == 0) return
            if (seqlen == 0) {
                status = "empty"
            } else {
                # getline returns -1 if the file cannot be opened, i.e. no output yet
                r = (getline tmp < outpath)
                close(outpath)
                status = (r >= 0) ? "done" : "pending"
            }
            print "q" n "\t" qseqid "\t" outpath "\t" status "\t" header > map
            if (status == "pending") {
                print ">q" n > pend
                printf "%s", seq > pend
            }
        }
        /^>/ {
            flush()
            n++
            header = $0
            seq = ""
            seqlen = 0

            # qseqid = first whitespace-delimited token, as BLAST would report it
            body = header
            sub(/^>[[:space:]]*/, "", body)
            split(body, w, /[[:space:]]+/)
            qseqid = (w[1] != "") ? w[1] : "unnamed_" n

            # Filename stem; disambiguate collisions after sanitising/truncation
            safe = sanitize(header)
            if (safe == "") safe = "unnamed_" n
            stem = safe
            k = 1
            while (safe in seen) { k++; safe = stem "_dup" k }
            seen[safe] = 1

            outpath = outdir "/" safe ".out"
            next
        }
        {
            line = $0
            gsub(/[-[:space:]]/, "", line)   # remove alignment gaps and stray whitespace/CR
            if (line != "") {
                seq = seq line "\n"
                seqlen += length(line)
            }
        }
        END {
            flush()
            printf "" > pend   # make sure the file exists even if nothing is pending
        }
    ' "$1"
}

# -----------------------------------------------------------------------------
# run_blast
#   Writes to OUT.part and renames it to OUT only if blastn exits successfully.
# -----------------------------------------------------------------------------
run_blast() {
    local query="$1" out="$2" threads="$3"
    rm -f "$out" "$out.part"
    # $BLAST_OPTIONS is deliberately unquoted so that it splits into separate arguments
    if blastn -query "$query" \
              -db "$DB_PATH" \
              -out "$out.part" \
              -outfmt "6 $BLAST_COLS" \
              -num_threads "$threads" \
              "${TAX_ARGS[@]}" \
              $BLAST_OPTIONS 2> "$out.err"; then
        mv -f "$out.part" "$out"
        return 0
    fi
    rm -f "$out.part"
    return 1
}

# -----------------------------------------------------------------------------
# demux
#   Splits one chunk's BLAST output into the per-sequence .out files. Synthetic
#   IDs in column 1 are replaced with the real qseqid. Members with no rows get
#   an empty .out (= no hits). Every file is written as .part and then renamed.
# -----------------------------------------------------------------------------
demux() {
    local mapf="$1" memf="$2" tabf="$3"
    awk -F'\t' -v OFS='\t' -v mapf="$mapf" -v memf="$memf" '
        FILENAME == mapf { qs[$1] = $2; op[$1] = $3; next }
        FILENAME == memf { if ($1 != "") { mem[$1] = 1; order[++nm] = $1 }; next }
        NF == 0 { next }
        !($1 in mem) {
            printf "demux: BLAST reported query id %s, which is not in this chunk\n", $1 > "/dev/stderr"
            bad = 1
            exit
        }
        { id = $1; $1 = qs[id]; rows[id] = rows[id] $0 "\n" }
        END {
            if (bad) exit 3
            for (i = 1; i <= nm; i++) {
                id = order[i]
                f = op[id] ".part"
                printf "%s", rows[id] > f
                close(f)
                print op[id]
            }
        }
    ' "$mapf" "$memf" "$tabf" > "$tabf.written" || return 1

    local p
    while IFS= read -r p; do
        mv -f "$p.part" "$p"
    done < "$tabf.written"
}

# -----------------------------------------------------------------------------
# run_chunk
#   Runs one chunk. If blastn or the demux fails, retries each member on its
#   own so that one bad record cannot lose the rest of the chunk. Sequences
#   that still fail get no .out file.
# -----------------------------------------------------------------------------
run_chunk() {
    local work="$1" k="$2" threads="$3"
    local query="$work/chunk_$k.fa"
    local members="$work/chunk_$k.members"
    local tab="$work/chunk_$k.tab"
    local nq t0=$SECONDS
    nq=$(wc -l < "$members")

    # Another chunk already hit a taxonomy error: don't start
    [ -e "$ABORT_FILE" ] && return 1

    if run_blast "$query" "$tab" "$threads"; then
        if demux "$work/map.tsv" "$members" "$tab"; then
            echo "  ✓ chunk $k ($nq sequences) done in $((SECONDS - t0))s"
            rm -f "$query" "$tab" "$tab.err" "$tab.written"
            return 0
        fi
    elif is_taxonomy_error "$tab.err"; then
        # Same error would happen for every sequence, so stop instead of retrying
        echo "  ✗ chunk $k: BLAST rejected the taxonomy filter: $(head -c 500 "$tab.err" | tr '\n' ' ')"
        cp "$tab.err" "$ABORT_FILE"
        rm -f "$query" "$tab" "$tab.err" "$tab.written"
        return 1
    fi

    echo "  ✗ chunk $k ($nq sequences) failed: $(head -c 500 "$tab.err" 2>/dev/null | tr '\n' ' ')"
    echo "    retrying chunk $k one sequence at a time..."

    local qid sq stab smem n_ok=0 n_fail=0
    while IFS= read -r qid; do
        [ -e "$ABORT_FILE" ] && break
        sq="$work/single_${qid}.fa"
        stab="$work/single_${qid}.tab"
        smem="$work/single_${qid}.members"
        awk -v id=">$qid" '/^>/ { p = ($0 == id) } p' "$query" > "$sq"
        echo "$qid" > "$smem"
        if run_blast "$sq" "$stab" "$threads" && demux "$work/map.tsv" "$smem" "$stab"; then
            ((n_ok++))
        else
            ((n_fail++))
            is_taxonomy_error "$stab.err" && cp "$stab.err" "$ABORT_FILE"
            echo "$qid" >> "$work/chunk_$k.failed"
            echo "    ✗ $qid failed on retry: $(head -c 300 "$stab.err" 2>/dev/null | tr '\n' ' ')"
        fi
        rm -f "$sq" "$smem" "$stab" "$stab.err" "$stab.written"
    done < "$members"

    echo "  chunk $k: recovered $n_ok/$nq sequences on retry ($n_fail still failing)"
    rm -f "$query" "$tab" "$tab.err" "$tab.written"
}

# -----------------------------------------------------------------------------
# finalise FASTA
#   Goes through the original input (gaps kept) and sends each record to
#   hits / no_hits / failed, based on its per-sequence .out file. Builds the
#   combined table in input order. Prints "n_hits n_no_hits n_failed".
# -----------------------------------------------------------------------------
finalise() {
    awk -v mapf="$2/map.tsv" -v hits="$3.part" -v nohits="$4.part" \
        -v comb="$5.part" -v failed="$6.part" -v cols="$BLAST_COLS" '
        BEGIN {
            while ((getline line < mapf) > 0) {
                split(line, a, "\t")
                qs[a[1]] = a[2]; op[a[1]] = a[3]; st[a[1]] = a[4]
            }
            close(mapf)
            h = cols
            gsub(/ /, "\t", h)
            print h > comb
            print "qseqid\treason\toriginal_header" > failed
            printf "" > hits
            printf "" > nohits
        }
        /^>/ {
            n++
            id = "q" n
            dest = ""
            if (st[id] == "empty") {
                reason = "no sequence left after removing gaps"
            } else {
                r = (getline row < op[id])
                if (r < 0) {
                    reason = "BLAST failed (no output)"
                } else if (r == 0) {
                    dest = nohits; nn++
                } else {
                    dest = hits; nh++
                    print row > comb
                    while ((getline row < op[id]) > 0) print row > comb
                }
                close(op[id])
            }
            if (dest == "") { nf++; print qs[id] "\t" reason "\t" $0 > failed }
            else print > dest
            next
        }
        dest != "" { print > dest }
        END { print nh + 0, nn + 0, nf + 0 }
    ' "$1"
}

# Running totals across input files
TOTAL_HITS=0
TOTAL_NO_HITS=0
TOTAL_FAILED=0

# -----------------------------------------------------------------------------
# process_file FASTA
# -----------------------------------------------------------------------------
process_file() {
    local fasta="$1"
    local filename base per_seq_dir work
    filename=$(basename "$fasta")
    base="${filename%.*}"
    per_seq_dir="${OUTPUT_DIR}/${base}"
    work="${WORK_ROOT}/${base}"
    mkdir -p "$per_seq_dir" "$work"

    # Leftovers from an interrupted run are never valid results
    find "$per_seq_dir" -maxdepth 1 -name '*.out.part' -delete

    local seq_count
    seq_count=$(grep -c '^>' "$fasta")
    if (( seq_count == 0 )); then
        echo "  ✗ Error: No FASTA sequences found in $fasta"
        return 1
    fi

    scan_fasta "$fasta" "$per_seq_dir" "$work"

    local n_done n_empty n_pending
    read -r n_done n_empty n_pending < <(awk -F'\t' '{ c[$4]++ } END { print c["done"] + 0, c["empty"] + 0, c["pending"] + 0 }' "$work/map.tsv")
    echo "  $seq_count sequences: $n_pending to BLAST, $n_done already done (resumed), $n_empty empty after removing gaps"

    if (( n_pending > 0 )); then
        # One chunk per CPU (fewer if there are fewer sequences); more chunks
        # only if MAX_CHUNK_SIZE requires it. Sequences are dealt out round-robin
        # so chunk sizes differ by at most one.
        local n_chunks min_chunks workers threads_per_blast
        n_chunks=$(( n_pending < THREADS ? n_pending : THREADS ))
        min_chunks=$(( (n_pending + MAX_CHUNK_SIZE - 1) / MAX_CHUNK_SIZE ))
        (( min_chunks > n_chunks )) && n_chunks=$min_chunks
        workers=$(( n_chunks < THREADS ? n_chunks : THREADS ))
        threads_per_blast=$(( THREADS / workers ))

        awk -v n="$n_chunks" -v dir="$work" '
            /^>/ {
                k = (c % n) + 1; c++
                f = dir "/chunk_" k ".fa"
                print > f
                print substr($0, 2) > (dir "/chunk_" k ".members")
                next
            }
            { print > f }
        ' "$work/pending.fa"

        echo "  Running $n_chunks chunk(s) of ~$(( (n_pending + n_chunks - 1) / n_chunks )) sequences: $workers concurrent blastn process(es) x $threads_per_blast thread(s)"

        local k
        for (( k = 1; k <= n_chunks; k++ )); do
            while (( $(jobs -rp | wc -l) >= workers )); do
                wait -n
            done
            [ -e "$ABORT_FILE" ] && break
            run_chunk "$work" "$k" "$threads_per_blast" &
        done
        wait

        if [ -e "$ABORT_FILE" ]; then
            echo "Error: BLAST rejected the taxonomy filter ($TAX_DESC). Stopping."
            echo "BLAST said: $(head -c 1000 "$ABORT_FILE" | tr '\n' ' ')"
            echo "Check the taxids exist in the database and that it is a v5 database with taxonomy."
            echo "Higher-level taxids need taxonomy4blast.sqlite3; otherwise use get_species_taxids.sh."
            exit 1
        fi
    fi

    local hit_file="${OUTPUT_DIR}/${base}_hits.fasta"
    local no_hit_file="${OUTPUT_DIR}/${base}_no_hits.fasta"
    local combined="${OUTPUT_DIR}/${base}_blast.tsv"
    local failed_file="${OUTPUT_DIR}/${base}_failed.tsv"
    local nh nn nf

    read -r nh nn nf < <(finalise "$fasta" "$work" "$hit_file" "$no_hit_file" "$combined" "$failed_file")
    mv -f "$hit_file.part" "$hit_file"
    mv -f "$no_hit_file.part" "$no_hit_file"
    mv -f "$combined.part" "$combined"
    mv -f "$failed_file.part" "$failed_file"

    if [ ! -s "$no_hit_file" ]; then
        echo "  ! No sequences passed the filter (no BLAST hits). Removing empty file."
        rm -f "$no_hit_file"
    else
        echo "  * $(basename "$no_hit_file"): $nn sequences without BLAST hits"
    fi

    if [ ! -s "$hit_file" ]; then
        echo "  ! No sequences had BLAST hits. Removing empty file."
        rm -f "$hit_file"
    else
        echo "  * $(basename "$hit_file"): $nh sequences with BLAST hits"
    fi
    echo "  * $(basename "$combined"): combined BLAST table"

    if (( nf > 0 )); then
        echo "  WARNING: $nf sequence(s) failed and are in NEITHER FASTA - see $(basename "$failed_file")"
    else
        rm -f "$failed_file"
    fi

    TOTAL_HITS=$(( TOTAL_HITS + nh ))
    TOTAL_NO_HITS=$(( TOTAL_NO_HITS + nn ))
    TOTAL_FAILED=$(( TOTAL_FAILED + nf ))
    rm -rf "$work"
}

# Check if input path exists
if [ ! -e "$INPUT_PATH" ]; then
    echo "Error: Input path $INPUT_PATH does not exist"
    exit 1
fi

# Display configuration
echo "===== BLAST Configuration ====="
echo "Input: $INPUT_PATH"
echo "Output directory: $OUTPUT_DIR"
echo "Database name: $DB_NAME"
echo "Database: $DB_PATH"
echo "BLAST options: $BLAST_OPTIONS"
echo "Taxonomy filter: $TAX_DESC"
echo "CPUs: $THREADS (max $MAX_CHUNK_SIZE sequences per blastn call)"
echo "Note: Dashes (-) are removed from sequences for BLAST searches, but preserved in both filtered outputs"
echo "============================="

if [ -d "$INPUT_PATH" ]; then
    echo "Input is a directory. Processing all FASTA files..."

    shopt -s nullglob
    fasta_files=( "$INPUT_PATH"/*.fa "$INPUT_PATH"/*.fasta "$INPUT_PATH"/*.fas )
    shopt -u nullglob

    if (( ${#fasta_files[@]} == 0 )); then
        echo "Error: No FASTA files found in $INPUT_PATH"
        exit 1
    fi
    echo "Found ${#fasta_files[@]} FASTA files to process"

    count=0
    for fasta_file in "${fasta_files[@]}"; do
        ((count++))
        echo "[$count/${#fasta_files[@]}] Processing file: $(basename "$fasta_file")"
        process_file "$fasta_file"
    done
else
    echo "Input is a single file: $(basename "$INPUT_PATH")"
    if [[ "$INPUT_PATH" != *.fa && "$INPUT_PATH" != *.fasta && "$INPUT_PATH" != *.fas ]]; then
        echo "Warning: Input file does not have a .fa/.fasta/.fas extension. Continuing anyway..."
    fi
    process_file "$INPUT_PATH" || exit 1
fi

echo "===== Results ====="
echo "Sequences with BLAST hits:    $TOTAL_HITS"
echo "Sequences without BLAST hits: $TOTAL_NO_HITS"
echo "Failed sequences:             $TOTAL_FAILED"
if (( TOTAL_FAILED > 0 )); then
    echo "WARNING: some sequences failed and are in neither FASTA (see *_failed.tsv)."
    echo "         Re-run the script to retry only those sequences."
fi
echo "Results saved to $OUTPUT_DIR"

# Calculate and print total runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(( (runtime % 3600) / 60 ))
seconds=$((runtime % 60))

echo "===== Runtime Summary ====="
echo "Total runtime: ${hours}h ${minutes}m ${seconds}s"
echo "Total seconds: ${runtime}"
echo "=========================="
