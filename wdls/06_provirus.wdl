version development

workflow MAG_provirus{

  input{
  ### Inputs to workflow 
  File vamboutput
  File checkvdb
  File Genecatalogue
  File TPMmatrix
  File workflowscripts
  File viral_results
  File MGVdb
  }

  ### Run CheckV in parallel 

  call divide_MAG_bins{
          input:
          vamboutput = vamboutput
  }
   
  Array[File] MAGbins = divide_MAG_bins.MAGbins
  Array[Int] binindex = range(length(MAGbins))

  scatter (i in binindex){
    
    File binlist = MAGbins[i]
    String part = i

  call checkv{
        input:
           vamboutput = vamboutput,
           binlist = binlist,
           checkvdb = checkvdb,
           workflowscripts = workflowscripts
    }
  }


  ### Identify most legit Viruses
  call select_viruses{
        input:
          checkv_viruses = checkv.checkv_viruses,
          checkv_completeness = checkv.checkv_completeness,
          checkv_quality = checkv.checkv_quality,
          workflowscripts = workflowscripts
  }

  call derepblast_viruses{
        input:
          viral_results = viral_results,
          checkv_HQ = select_viruses.checkv_HQ,
          Genecatalogue = Genecatalogue,
          workflowscripts = workflowscripts
  }

  call getvirusTPM{
      input:
        virusproducts= derepblast_viruses.virusproducts,
        TPMmatrix = TPMmatrix,
        workflowscripts = workflowscripts
  }

  call virus_MGV_clustering{
      input:
        virusproducts= derepblast_viruses.virusproducts,
        workflowscripts = workflowscripts,
        MGVdb = MGVdb
  }


  output {
    File Provirusresults = getvirusTPM.ProVirusResults
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

task checkv{
  input{
    File vamboutput
    File binlist
    File checkvdb
    File workflowscripts
  }
  command{
      tar -xvf ${workflowscripts}
      tar -xvf ${vamboutput}
      mkdir -p MAGs 
      for file in `cat ~{binlist}`; do 
          mv $file MAGs
      done
      cat MAGs/*fna.gz > all_MAG_contigs.fna.gz 

      ### Run checkv
      tar -xvf ${checkvdb}
      export CHECKVDB=checkv-db-v1.0

      pigz -d all_MAG_contigs.fna.gz
      python pyscripts/fastx_rename.py --i all_MAG_contigs.fna --o all_MAG_contigs.renamed.fna

      checkv end_to_end all_MAG_contigs.renamed.fna checkvrun -t 34

      ### Export relevant files
      mv checkvrun/completeness.tsv completeness.MAG.tsv
      mv checkvrun/quality_summary.tsv quality_summary.MAG.tsv
      mv checkvrun/viruses.fna viruses.fna
      pigz viruses.fna
  }
  output{
    File checkv_viruses = 'viruses.fna.gz'
    File checkv_completeness = 'completeness.MAG.tsv'
    File checkv_quality = 'quality_summary.MAG.tsv'
  }
  runtime{
      docker: "quay.io/joacjo/vamb_virus"
      disks: "local-disk 250 HDD"
      memory: "30GB"
      cpu : "34"
  }
}

task select_viruses{
   input{
    Array[File] checkv_viruses
    Array[File] checkv_completeness
    Array[File] checkv_quality
    File workflowscripts
  }
  command <<<
      tar -xvf ~{workflowscripts}

      cat ~{sep=' ' checkv_viruses}  > proviruses.fna.gz
      cat ~{sep=' ' checkv_completeness}  > completeness.MAG.tsv
      cat ~{sep=' ' checkv_quality} > quality_summary.MAG.tsv

      Rscript --vanilla --slave Rscripts/parse_checkv_results.R quality_summary.MAG.tsv completeness.MAG.tsv proviruses
      cat proviruses.HQ.tsv proviruses.MQ.tsv proviruses.LQ.tsv proviruses.provirus.tsv > HQMQLQ.proviruses.txt

      ### Pack results 
      tar -czvf checkv_selected_viruses.tar.gz proviruses.fna.gz HQMQLQ.proviruses.txt quality_summary.MAG.tsv completeness.MAG.tsv
  >>>
  output{
    File checkv_HQ = "checkv_selected_viruses.tar.gz"
  }
  runtime{
      docker: "gcr.io/microbiome-xavier/r-docker:010319"
      disks: "local-disk 250 HDD"
      cpu : "1"
      preemptible: 2
  }
}


task derepblast_viruses{
  input{
    File viral_results  
    File checkv_HQ
    File Genecatalogue
    File workflowscripts
  }
  command{
      ### 
      tar -xvf ${checkv_HQ}
      tar -xvf ${viral_results}
      cut -f1 HQMQLQ.proviruses.txt > HQMQLQ.names
      seqtk subseq proviruses.fna.gz HQMQLQ.names > HQMQLQ.proviruses.fna

      ### Concat MAG-derived proviruses with Viral-bins
      cat virusresults/viralreps.fna HQMQLQ.proviruses.fna > HQMQLQ.allviruses.fna 

      ### 
      makeblastdb -in HQMQLQ.allviruses.fna -dbtype nucl -out HQMQLQdb
      blastn -query HQMQLQ.allviruses.fna -db HQMQLQdb -outfmt '6 std qlen slen' \
      -max_target_seqs 10000 -perc_identity 90 \
      -out virusxvirus.m6 -num_threads 40

      python /app/anicalc.py -i virusxvirus.m6 -o virusxvirus.ani.tsv
      python /app/aniclust.py --fna HQMQLQ.allviruses.fna --ani virusxvirus.ani.tsv --out virus.clusters --min_ani 95 --min_qcov 0 --min_tcov 85
      
      cut -f1 virus.clusters  > viralreps.txt
      seqtk subseq HQMQLQ.allviruses.fna viralreps.txt > viralreps.allviruses.fna 
      makeblastdb -in viralreps.allviruses.fna -dbtype nucl -out HQMQLQdbreps

      ### Blast Genecatalogue against the viral references
      gunzip -dc ${Genecatalogue} > MGX.genecat.fna

      blastn -task megablast -db HQMQLQdbreps -query MGX.genecat.fna -out HQMQLQ.allviruses.genecat.m6 -outfmt 6 -num_threads 40 \
      -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1

      tar -czvf provirus.clustering.tar.gz HQMQLQ.proviruses.txt viralreps.allviruses.fna HQMQLQ.proviruses.fna virus.clusters HQMQLQ.allviruses.genecat.m6
  }
  output{
    File virusproducts = "provirus.clustering.tar.gz"
  }
  runtime{
      docker: "quay.io/joacjo/vamb_virus:latest"
      disks: "local-disk 250 HDD"
      memory: "60GB"
      cpu : "40"
      preemptible: 1
  }
}

task getvirusTPM{
    input{
        File virusproducts
        File TPMmatrix
        File workflowscripts
    }
    command{
        ### 
        tar -xvf ${virusproducts}
        tar -xvf ${workflowscripts}
        mkdir -p matrices
        gunzip -fdc ${TPMmatrix} > TPMmatrix.txt

        install2.r --error fastmatch plyr
        Rscript --vanilla --slave Rscripts/get_viral_abundance_matrix.R HQMQLQ.allviruses.genecat.m6 TPMmatrix.txt matrices

        ### pack results
        mkdir ProVirusresults
        mv matrices/* ProVirusresults
        mv viralreps.allviruses.fna ProVirusresults 
        mv virus.clusters ProVirusresults
        mv HQMQLQ.allviruses.genecat.m6 ProVirusresults
        tar -czvf Provirus.matrices.tar.gz ProVirusresults
    }
    output{
    File ProVirusResults = "Provirus.matrices.tar.gz"
    }
    runtime{
      docker: "rocker/tidyverse:latest"
      disks: "local-disk 250 HDD"
      memory: "60GB"
      cpu : "4"
      preemptible: 1
    }
}


task virus_MGV_clustering{
    input{
      File virusproducts
      File workflowscripts
      File MGVdb
    }
    command <<<
      tar -xvf ~{virusproducts}
      tar -xvf ~{workflowscripts}
      
      ### Combine with MGVdb 
      pigz -dc ~{MGVdb} > mgvreps.fna
      cat viralreps.allviruses.fna mgvreps.fna > viruses.fna

      makeblastdb -in viruses.fna -dbtype nucl -out virusesDB
      blastn -query viruses.fna -db virusesDB -outfmt '6 std qlen slen' \
      -max_target_seqs 10000 -perc_identity 90 \
      -out virusxvirus.m6 -num_threads 40

      python /app/anicalc.py -i virusxvirus.m6 -o virusxvirus.ani.tsv
      python /app/aniclust.py --fna viruses.fna --ani virusxvirus.ani.tsv --out virus.clusters.95 --min_ani 95 --min_qcov 0 --min_tcov 70
      python /app/aniclust.py --fna viruses.fna --ani virusxvirus.ani.tsv --out virus.clusters.80 --min_ani 80 --min_qcov 0 --min_tcov 70

      mkdir MGVclustering 
      mv virusxvirus.ani.tsv MGVclustering
      mv virus.clusters.95 MGVclustering
      mv virus.clusters.80 MGVclustering
      mv virusxvirus.m6 MGVclustering
      pigz MGVclustering/virusxvirus.m6
      tar -czvf ProvirusMGVclustering.tar.gz MGVclustering
    >>>
    output{
    File MGVcluster = "ProvirusMGVclustering.tar.gz"
    }
    runtime{
    docker: "quay.io/joacjo/vamb_virus"
    disks: "local-disk 250 HDD"
    memory: "80GB"
    cpu : "40"
    preemptible: 1
  }

}
