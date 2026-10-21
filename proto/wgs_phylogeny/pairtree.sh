
####
    BaseDir="/data/wgs/RUO_Service"
    BatchID="WGS-26-02"
    AnalysisList=${BaseDir}/${BatchID}/meta/clonality.evolution.analysis.list
####
####
    PairtreeDir="${BaseDir}/${BatchID}/analysis/cancer-evolution-pairtree2"
    mkdir -p "${PairtreeDir}"
    cd ${PairtreeDir}
####

    singularity exec -B /storage,/data /storage/images/pairtree.sif python3 /opt/conda/pairtree/util/fix_bad_var_read_prob.py \
        ${PairtreeDir}/${BatchID}.pairtree.ssm.txt \
        ${PairtreeDir}/${BatchID}.pairtree.params.json \
        ${PairtreeDir}/${BatchID}.pairtree.ssm.fixed.txt \
        ${PairtreeDir}/${BatchID}.pairtree.params.fixed.json

    singularity exec -B /storage,/data /storage/images/pairtree.sif python3 /opt/conda/pairtree/bin/clustervars \
        ${PairtreeDir}/${BatchID}.pairtree.ssm.fixed.txt \
        ${PairtreeDir}/${BatchID}.pairtree.params.fixed.json \
        ${PairtreeDir}/${BatchID}.pairtree.params.fixed.json

    singularity exec -B /storage,/data /storage/images/pairtree.sif python3 /opt/conda/pairtree/bin/pairtree \
        --params ${PairtreeDir}/${BatchID}.pairtree.params.fixed.json \
        --parallel 24 \
        --seed=5555 \
        ${PairtreeDir}/${BatchID}.pairtree.ssm.fixed.txt \
        ${PairtreeDir}/${BatchID}.pairtree.result.npz

    singularity exec -B /storage,/data /storage/images/pairtree.sif python3 /opt/conda/pairtree/bin/plottree \
        ${PairtreeDir}/${BatchID}.pairtree.ssm.fixed.txt \
        ${PairtreeDir}/${BatchID}.pairtree.params.fixed.json \
        ${PairtreeDir}/${BatchID}.pairtree.result.npz \
        ${PairtreeDir}/${BatchID}.pairtree.result.plots.html

    python3 /storage/home/kangsm/apps/pairtree/extract_pairtree_node_variants.py ${PairtreeDir}/${BatchID}.pairtree.result.npz ${PairtreeDir}/${BatchID}.best.tree.node.variants.tsv
