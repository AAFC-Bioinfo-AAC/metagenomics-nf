/* 
This workflow renames sequences by sample id given a map_file
*/

include { RENAME_SEQUENCES } from '../../modules/local/rename_sequences'

workflow RENAME {
    
    println "You are using the *RENAME* subworkflow."

    take: 
        data
        map_file

    main:
      RENAME_SEQUENCES(data.collect(), map_file)
      
    emit:
      RENAME_SEQUENCES.out

}