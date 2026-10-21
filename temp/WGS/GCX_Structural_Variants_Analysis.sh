





    RunSingularity="singularity exec -B /storage,/data /storage/images"
    ReferenceFasta=/storage/references_and_index/hg38/fasta/cbNIPT/hg38.fa
    Threads=
### breakseq2

    ${RunSingularity}/breakseq2-2.2.sif run_breakseq2.py \
        --bams
        --nthreads ${Threads}
        --chromosomes
        --work
        --reference
        --bplib_gff 
        --min_span
        --window
        --min_overlap
        --junction_length


### breakdancer
    config 
    readgroup:<RGID> \t platform:illumina \t map:<aa.bam> \t readlen: \t lib: \t num: \t lower: \t upper: \t mean: \t SWnormality: \t exe:
    
    ${RunSingularity}/breakdancer-1.4.5.sif bam2cfg \
        --input-file
        --min
        --min-mapq
        --output
        --dist-output


    ${RunSingularity}/breakdancer-1.4.5.sif breakdancer-max \
        config 
        -o chr
        -s 5 \
        -q 10
        -r 2
        -a -h



    ${RunSingularity}/delly_v1.7.2.sif delly call --help
        --svtype INV 
        --genome ${ReferenceFasta}
        --outfile 
        --threads ${Threads}
        --map-qual 10
        --qual-tra 
        --minclip 
        --minrefsep
        --maxreadsep
        --max-reads
        --geno-qual


