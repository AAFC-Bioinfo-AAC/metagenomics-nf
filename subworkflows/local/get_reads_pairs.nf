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
       .filter( ~/^.*fastq.gz/ )
       .map { file ->
          def key_part1 = file.name.toString().tokenize('_').get(0)
          key = key_part1
          return tuple(key, file)
       }
       .groupTuple(size:2)
       .flatten()
       .collate ( 3 )
       .set { read_pairs_ch }
      
    emit:
      read_pairs_ch

}