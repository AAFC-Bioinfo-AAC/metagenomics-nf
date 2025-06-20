/*
This workflow creates a `folder_ch` channel that emits tuples containing 2 
elements: the sample ID, and the path containing a folder with stuff.
*/

workflow GET_FOLDER {
    
    println "You are using the *GET_FOLDER* subworkflow."

    take: data
    
    main:
       data.flatten()
       .map { dir ->
          def key_part1 = dir.name.toString()
          key = key_part1
          return tuple(key, dir)
       }
       //.view()       
       .set { folder_ch }

    emit:
      folder_ch
      
}