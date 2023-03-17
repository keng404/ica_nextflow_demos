def VERSION = '1.0' // Version information not provided by tool on CLI
process GET_MLST_SRST2 {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/python:2.7--1"
    input:
    tuple val(meta),  path(taxonomy), val(status)
    output:
    tuple val(meta), path("*_getMLST_out.txt")                                 , optional:true, emit: getMLSTs
    tuple val(meta), path("*.fasta")                                           , optional:true, emit: fastas
    tuple val(meta), path("*_profiles.csv")                                    , optional:true, emit: profiles
    tuple val(meta), path("*_getMLST_out_temp.txt")                            , optional:true, emit: getMLSTs_checker
    tuple val(meta), path("*_temp.fasta")                                      , optional:true, emit: fastas_checker
    tuple val(meta), path("*_profiles_temp.csv")                               , optional:true, emit: profiles_checker
    tuple val(meta), path("*_pull_dates.txt")                                  , optional:true, emit: pull_date
    path "versions.yml"                                                                       , emit: versions
    when:
    (task.ext.when == null || task.ext.when) //& "${status[0]}" == "False"
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${status[0]}" == "False" ]]; then
      genus="empty"
      species="empty"
      today=\$(date '+%Y-%m-%d')
      test_title=\$(tail -n2 ${taxonomy} | head -n1)
      echo "-\${test_title}-"
      if [[ \$test_title = "G:"* ]]; then
        species=\$(tail -n1 ${taxonomy} | cut -f2)
        genus=\$(tail -n2 ${taxonomy} | head -n1 | cut -f2)
      elif [[ \$test_title = "s:"* ]]; then
        species=\$(tail -n2 ${taxonomy} | head -n1 | cut -f2)
        genus=\$(tail -n3 ${taxonomy} | head -n1 | cut -f2)
      else
        echo "-\$test_title-"
      fi
      echo "\${genus}___\${species}"
      python ${workflow.launchDir}/bin/convert_taxonomy_with_complexes_to_pubMLST.py --genus "\${genus}" --species "\${species}" > DB_defs.txt
      dbline=\$(tail -n1 DB_defs.txt)
      echo "\$dbline"
      IFS=',' read -r -a db_array <<< "\$dbline"
      echo "\${#db_array[@]}-\${db_array[@]}"
      #num_dbs="\${#db_array[@]}"
      counter=1
      for entry in "\${db_array[@]}"
      do
        echo "Entry#\${counter}-\${entry}|"
        entry_no_spaces="\${entry// /_}"
        if [[ "\${entry}" = "No match found" ]]; then
          touch "\${entry_no_spaces}.fasta"
          touch "\${entry_no_spaces}_profiles.csv"
          touch "\${entry_no_spaces}_pull_dates.txt"
          touch "\${entry_no_spaces}_temp.fasta"
          touch "\${entry_no_spaces}_profiles_temp.csv"
          echo "DB:No match found(\${genus} \${species})       defs:\${entry_no_spaces}_profiles.csv        del:''" > \${entry_no_spaces}_getMLST_out.txt
          cp "\${entry_no_spaces}_getMLST_out.txt" "\${entry_no_spaces}_getMLST_out_temp.txt"
        else
          if [[ "\${entry}" = "Streptococcus thermophilus" ]]; then
            python ${workflow.launchDir}/bin/getMLST2_phoenix.py --species "\$entry" --force_scheme_name
          else
            python ${workflow.launchDir}/bin/getMLST2_phoenix.py --species "\$entry"
          fi
          if [[ ! -f dbases.xml ]]; then
            touch "\${entry_no_spaces}.fasta"
            touch "\${entry_no_spaces}_profiles.csv"
            touch "\${entry_no_spaces}_temp.fasta"
            touch "\${entry_no_spaces}_profiles_temp.csv"
            echo "DB:Server down(\${genus} \${species})       defs:\${entry_no_spaces}_profiles.csv        del:''" > \${entry_no_spaces}_getMLST_out.txt
            cp "\${entry_no_spaces}_getMLST_out.txt" "\${entry_no_spaces}_getMLST_out_temp.txt"
          else
            if [[ "\${entry}" = *"baumannii#1" ]]; then
              sed -i -e 's/Oxf_//g' "\${entry_no_spaces}.fasta"
              sed -i -e 's/Oxf_//g' "\${entry_no_spaces}_profiles.csv"
            elif [[ "\${entry}" = *"baumannii#2" ]]; then
              sed -i -e 's/Pas_//g' "\${entry_no_spaces}.fasta"
              sed -i -e 's/Pas_//g' "\${entry_no_spaces}_profiles.csv"
            fi
            cp "\${entry_no_spaces}.fasta" "\${entry_no_spaces}_temp.fasta"
            cp "\${entry_no_spaces}_profiles.csv" "\${entry_no_spaces}_profiles_temp.csv"
            cp "\${entry_no_spaces}_getMLST_out.txt" "\${entry_no_spaces}_getMLST_out_temp.txt"
          fi
          echo "\${today}" > "\${entry_no_spaces}_pull_dates.txt"
        fi
        counter=\$(( counter + 1))
      done
    else
        echo "Did not run" > "${prefix}_getMLST_out_temp.txt"
        echo "Did not run" > "${prefix}_temp.fasta"
        echo "Did not run" > "${prefix}_profiles_temp.csv"
    fi
    echo -e "\"${task.process}\":
        getMLST: $VERSION" > versions.yml
    """
}
