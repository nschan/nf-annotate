/*
 Below is a somewhat parallel approach to running augustus
 This is far from optimal, it just shoots all subjobs into the scheduler at the same time,
 this may or may not work so well.
 I have not managed to implement something that runs this batched, but most of the jobs should be fast?
 It is largely based on https://bioinf.uni-greifswald.de/bioinf/wiki/pmwiki.php?n=Augustus.ParallelPred
 */

process AUGUSTUS_PARALLEL {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(accession_genome)
        val(species)

    output:
        tuple val(meta), path("*_augustus.gff"), emit: augustus_gff

    script:
        """
        splitMfasta.pl ${accession_genome}
        for f in *.split.*; do
            NAME=`grep ">" \$f | sed 's/>//'`; mv \$f \${NAME}.fa &
        done
        wait
        summarizeACGTcontent.pl ${accession_genome} > summary.out
        grep "bases" summary.out | awk '{print \$3".fa","1",\$1}' > chr.txt
        createAugustusJoblist.pl \\
            --sequences=chr.txt \\
            --wrap="#" \\
            --overlap=5000 \\
            --chunksize=1252500 \\
            --outputdir="${meta}" \\
            --joblist=jobs.lst \\
            --jobprefix="" \\
            --command "augustus --exonnames=on --codingseq=on --species=${species}"
        nJobs="\$(tail -n1 jobs.lst)"
        mkdir ${meta}
        for ((j=nJobs; j>0; j--)); do 
                cat \$j | bash & 
        done 
        wait
        cat ${meta}/*gff | join_aug_pred.pl > ${meta}_augustus.gff.tmp
        awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "AUGUSTUS"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${meta}_augustus.gff.tmp > ${meta}_augustus.gff
        rm ${meta}_augustus.gff.tmp
        rm -r ${meta}
        cat jobs.lst | xargs rm
        """
}

/*
Here is the classic single non-parallel approach. 
It is much slower than the module above...
*/

process AUGUSTUS {
    tag "$meta"
    label 'process_medium'
    publishDir(
      path: { "${params.out}/${task.process}".replace(':','/').toLowerCase() }, 
      mode: 'copy',
      overwrite: true,
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    ) 
    input:
        tuple val(meta), path(accession_genome)
        val(species)

    output:
        tuple val(meta), path("*_augustus.gff"), emit: augustus_gff

    script:
        """
        augustus \\
        --species=${species} \\
        ${accession_genome} > ${meta}_augustus.gff.tmp
        awk 'BEGIN {OFS="\\t"}; /^[^#]/ {\$2 = "AUGUSTUS"; nine=\$9
             for(i=10; i <= NF; i++) nine=nine" "\$i
             print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, nine}' ${meta}_augustus.gff.tmp > ${meta}_augustus.gff
        rm ${meta}_augustus.gff.tmp
        """
}
 
 
 