version development

workflow MAG_quality_tax{

  input{
  ### Inputs to workflow 
  File vamboutput
  File workflowscripts
  File checkm_db
  File gtdbtk_db
  File CRTJAR

  }

  call divide_MAG_bins{
          input:
          vamboutput = vamboutput
  }
   
  Array[File] MAGbins = divide_MAG_bins.MAGbins
  Array[Int] binindex = range(length(MAGbins))

  scatter (i in binindex){
    
    File binlist = MAGbins[i]
    String part = i
  
    call checkm{
        input:
           vamboutput = vamboutput,
           binlist = binlist,
           part = part,
           checkmdb = checkm_db
    }

    call cctyper{
        input:
          vamboutput = vamboutput,
          binlist = binlist,
          workflowscripts = workflowscripts,
          CRTJAR = CRTJAR,
          part = part

    }

    call gtdbtk{
        input:
           vamboutput = vamboutput,
           MQNCMAGs = checkm.MQNCMAGs,
           part = part,
           gtdbtk_db = gtdbtk_db
    }
  }

  call collect_results{
        input:
          checkmfiles = checkm.checkm_file,
          spacers = cctyper.MAGspacers,
          crisprs = cctyper.MAGcrisprs,
          gtdbtk = gtdbtk.gtdbtk,
          crtspacers = cctyper.MAGcrtspacers
          
  }

  call format_taxonomy_table{
        input:
        gtdbktax = collect_results.gtdbktax
        workflowscripts = workflowscripts
  }

  output {
    File checkm_combined = collect_results.checkm_combined
    File MAGcrisprs = collect_results.MAGcrisprs
    File Gtdbtktax = collect_results.gtdbktax
    File Taxonomy = format_taxonomy_table.taxformatted
  }
}

### Write concatenated VAMB bins (of contigs)
task divide_MAG_bins{
  input{
    File vamboutput
  }
  command{
    tar -xvf ${vamboutput}

    ### Assuming that bins are splitted already by sample and located under vamb/bins
    ### We can just divide the directory of bins into smaller chunks 

    ls vamb/bins/*fna.gz > all_MAG_bins.txt
    mkdir binfiles

    ### SPlit into 500 MAG bin chunks that can be parsed seperately
    split -l 500 all_MAG_bins.txt binfiles/split

  }
  output{
    Array[File] MAGbins = glob("binfiles/split*")
  }

  runtime{
      docker: "ubuntu:16.04"
      disks: "local-disk 150 HDD"
      cpu : "1"
      preemptible: 1
  }
}

### Only bacterial and eukaryotic completeness estimates
task checkm{
  input{
    File vamboutput
    File binlist
    String part
    File checkmdb
  }
  command <<<
    tar -xvf ~{vamboutput}
    mkdir -p MAGs 
    for file in `cat ~{binlist}`; do 
        mv $file MAGs
    done

    ### CheckM database
    mkdir checkm_database
    tar -xvf ~{checkmdb} -C checkm_database
    checkm data setRoot checkm_database

    ### RUn CheckM
    mkdir -p checkmout
    pigz -d MAGs/*gz
    checkm lineage_wf -t 40 ./MAGs ./checkmout > checkmout/checkm.output.~{part} 
    grep "::" checkmout/checkm.output.~{part}  > checkmout/checkm.file
    
    awk '$13>=90 && $14<=5' checkmout/checkm.file > checkmout/NC.checkm.output.~{part} 
    awk '$13>=50 && $14<=10' checkmout/checkm.file > checkmout/MQNC.checkm.output.~{part}
  >>>
  output{
    File NCMAGs = "checkmout/NC.checkm.output.~{part}"
    File MQNCMAGs = "checkmout/MQNC.checkm.output.~{part}"
    File checkm_file = "checkmout/checkm.output.~{part}"
  }
  runtime{
      docker: "biojokke/checkm:1.1.3"
      disks: "local-disk 400 HDD"
      memory: "80GB"
      cpu : "40"
      preemptible: 1
  }
}

task gtdbtk{
  input{
    File vamboutput
    File MQNCMAGs
    String part
    File gtdbtk_db
  }
  command <<<
    tar -xvf ~{vamboutput}
    mkdir -p MAGs gtdbtk scratch results

    for MAG in `awk '{print $1}' ~{MQNCMAGs}`; do 
        cp vamb/bins/${MAG}.fna.gz MAGs
    done
    
    tar -xvf ~{gtdbtk_db}
    export GTDBTK_DATA_PATH=release202
    gtdbtk classify_wf \
        --scratch_dir scratch \
        --pplacer_cpus 1 \
        --extension gz \
        --genome_dir MAGs \
        --out_dir gtdbtk \
        --cpus 40 \
        --prefix gtdbtk \
        --force
    ls gtdbtk
    cp gtdbtk/gtdbtk.bac120.markers_summary.tsv results/~{part}.bac120.markers_summary.tsv
    cp gtdbtk/gtdbtk.ar122.markers_summary.tsv results/~{part}.ar122.markers_summary.tsv
    cp gtdbtk/gtdbtk.bac120.markers_summary.tsv results/~{part}.bac120.markers_summary.tsv
    cp gtdbtk/gtdbtk.translation_table_summary.tsv results/~{part}.translation_table_summary.tsv
    cp gtdbtk/gtdbtk.bac120.classify.tree results/~{part}.bac120.classify.tree
    cp gtdbtk/gtdbtk.ar122.summary.tsv results/~{part}.gtdbtk.ar122.summary.tsv
    cp gtdbtk/gtdbtk.bac120.summary.tsv results/~{part}.gtdbtk.bac120.summary.tsv
    tar -czvf gtdbtk.~{part}.tar.gz results
  >>>
  output{
    File gtdbtk = "gtdbtk.~{part}.tar.gz"
  }
  runtime{
      docker: "ecogenomic/gtdbtk"
      disks: "local-disk 500 HDD"
      memory: "160GB"
      cpu : "40"
      preemptible: 1
  }
}

### CCtyper - CRISPRs in MAGss
### Also run CRT
task cctyper{
  input{
    File vamboutput
    File binlist
    File workflowscripts
    File CRTJAR
    String part
  }
  command <<<
    tar -xvf ~{workflowscripts}
    tar -xvf ~{vamboutput}
    mkdir -p MAGs spacers predictions mkdir crtspacers crt
    export CCTYPER_DB="/miniconda/db"

    for MAG in `awk '{print $1}' ~{binlist}`; do 
        gunzip vamb/bins/${MAG}.fna.gz
        python pyscripts/rename_contigs.py --i vamb/bins/${MAG}.fna --o vamb/bins/${MAG}.renamed.fna --id ${MAG}
        cctyper vamb/bins/${MAG}.renamed.fna tmp --prodigal meta -t 20 --db /miniconda/db
        cat tmp/spacers/* > tmp/all.spacers.fna
        mv tmp/all.spacers.fna spacers/${MAG}.spacers.fna
        mv tmp/crisprs_all.tab predictions/${MAG}.tab
        rm -rf tmp

        java -cp ~{CRTJAR} crt -minNR 2 vamb/bins/${MAG}.fna ${MAG}.crt
        python pyscripts/CRT_extract.py -i ${MAG}.crt -o ${MAG}.crtspacers.fna
        mv ${MAG}.crtspacers.fna crtspacers/${MAG}.crtspacers.fna
        mv ${MAG}.crt crt
    done
    cat crtspacers/* > MAG.crisprs.crt.~{part}.fna
    cat predictions/* > MAG.crisprs.~{part}.tab
    cat spacers/* > MAG.spacers.~{part}.fna

  >>>
  output{
    File MAGspacers = "MAG.spacers.~{part}.fna"
    File MAGcrisprs = "MAG.crisprs.~{part}.tab"
    File MAGcrtspacers = "MAG.crisprs.crt.~{part}.fna"
  }
  runtime{
      docker: "quay.io/joacjo/cctyper:latest"
      disks: "local-disk 200 HDD"
      memory: "40GB"
      cpu: "20"
      preemptible: 1
  }
}


task collect_results{
  input{
    Array[File] checkmfiles
    Array[File] spacers
    Array[File] crisprs
    Array[File] gtdbtk
    Array[File] crtspacers

  }
  command <<<
    cat ~{sep=' ' crtspacers} > all.crtspacers.fna
    cat ~{sep=' ' checkmfiles} > all.checkm.tsv
    cat ~{sep=' ' spacers} > all.spacers.fna
    cat ~{sep=' ' crisprs} > all.crisprs.tab
    tar -czvf MAGcrisprs.tar.gz all.spacers.fna all.crisprs.tab all.crtspacers.fna

    ### Unpack GTDBTK files
    mkdir taxonomyfiles
    for f in ~{sep=' ' gtdbtk}; do
      tar -xvf ${f}
      mv results/*.summary.tsv taxonomyfiles
      rm -rf results
    done
    cat taxonomyfiles/* > gtdbtk.taxonommy.tsv

  >>>
  output{
    File checkm_combined = "all.checkm.tsv"
    File MAGcrisprs = "MAGcrisprs.tar.gz"
    File gtdbktax = "gtdbtk.taxonommy.tsv"
  }
  runtime{
      docker: "ubuntu:16.04"
      disks: "local-disk 150 HDD"
      cpu : "1"
      preemptible: 1
  }

}

task format_taxonomy_table{
  input{
    File gtdbktax
    File workflowscripts
  }
    command <<<
    python pyscripts/format_gtdbtktax.py -m gtdbtk.taxonommy.tsv -o gtdbtax.formatted.txt 
    >>>
  output{
    File taxformatted = "gtdbtax.formatted.txt"
  }
  runtime{
      docker: "quay.io/joacjo/cctyper:latest"
      disks: "local-disk 150 HDD"
      cpu : "1"
      preemptible: 1
  }
} 