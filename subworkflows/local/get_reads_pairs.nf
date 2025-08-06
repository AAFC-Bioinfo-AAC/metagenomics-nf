/*
This workflow creates the `read_pairs_ch` channel that emits tuples containing
three elements: the pair ID, the first read-pair file and the second read-pair
file.
*/

workflow GET_READS_PAIRS {
    
    println "You are using the *GET_READS_PAIRS* subworkflow."

    take: data
    
    main:
       data.flatten()
       .filter( ~/^.*R1.fastq.gz/ )
       .map { file ->
          def key = file.name.toString().tokenize('_').get(0)
          return tuple(key, file)
       }
       .set { R1_ch }
      
       data.flatten()
       .filter( ~/^.*R2.fastq.gz/ )
       .map { file ->
          def key = file.name.toString().tokenize('_').get(0)
          return tuple(key, file)
       }
       .set { R2_ch }

      R1_ch.join(R2_ch)
      .set { read_pairs_ch}

    emit:
      read_pairs_ch

}