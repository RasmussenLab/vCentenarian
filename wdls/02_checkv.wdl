version development

workflow checkv{

  input{
  ### Inputs to workflow 
  File allcontigs
  File allcontigs_fai
  File vamboutput
  File checkvdb
  File workflowscripts
  String vamb_split

  }

  call write_smallbins{
          input:
          allcontigs = allcontigs,
          vamboutput = vamboutput,
          vamb_split = vamb_split,
          workflowscripts = workflowscripts
  }
   
  Array[File] binarray = write_smallbins.smallbins
  Array[Int] binindex = range(length(binarray))

  scatter (i in binindex){
    
    File binfile = binarray[i]
    String part = i
  
    call checkv{
        input:
           binfile = binfile,
           checkvdb = checkvdb,
           part = part
    }
  }
  
  ### Organise bins by quality based on CheckV results
  call select_viruses{
        input:
          quality_summaries = checkv.checkv_qual,
          completeness = checkv.checkv_comp,
          workflowscripts = workflowscripts
  }
  call pack_viruses{
        input:
          viralbins = checkv.checkv_fasta,
          checkv_HQ = select_viruses.checkv_HQ,
          checkv_MQ = select_viruses.checkv_MQ,
          checkv_LQ = select_viruses.checkv_LQ
  }

  output {
    File CheckvFiles = pack_viruses.Checkvresults
  }
}


### Write concatenated VAMB bins (of contigs)
task write_smallbins{
  input{
    File allcontigs
    File vamboutput
    String vamb_split
    File workflowscripts
  }
  command{
    tar -xvf ${vamboutput}
    tar -xvf ${workflowscripts}

    mkdir -p smallbins
    python /app/write_concat_bins.py -f ${allcontigs} -o smallbins -c vamb/clusters.tsv -s ${vamb_split}
    
    ### 
    gzip smallbins/*fna
    tar -zcvf smallvamb.tar.gz smallbins

  }
  output{
    Array[File] smallbins = glob("smallbins/*fna.gz")
    File smallbins_tar = "smallvamb.tar.gz" 
  }

  runtime{
      docker: "quay.io/joacjo/vamb"
      disks: "local-disk 150 HDD"
      memory: "80GB"
      cpu : "1"
  }
}

task checkv{
  input{
    File binfile
    File checkvdb
    String part
    
  }
  command{
      ### Run checkv
      tar -xvf ${checkvdb}
      export CHECKVDB=checkv-db-v1.0
      pigz -dc ${binfile}  > binfile.fna
      checkv end_to_end binfile.fna checkvrun -t 24

      ### Export relevant files
      mv checkvrun/completeness.tsv completeness.${part}.tsv
      mv checkvrun/quality_summary.tsv quality_summary.${part}.tsv
      mv checkvrun/contamination.tsv contamination.${part}.tsv
      mv checkvrun/viruses.fna viruses.${part}.fna
      pigz viruses.${part}.fna

  }
  output{
    File checkv_comp = "completeness.${part}.tsv"
    File checkv_qual = "quality_summary.${part}.tsv"
    File checkv_cont = "contamination.${part}.tsv"
    File checkv_fasta = "viruses.${part}.fna.gz"
  }

  runtime{
      docker: "quay.io/joacjo/vamb_virus"
      disks: "local-disk 250 HDD"
      memory: "40GB"
      cpu : "24"
      preemptible: 1
  }
}

task select_viruses{
   input{
    Array[File] quality_summaries
    Array[File] completeness 
    File workflowscripts
  }
  command{
      cat ${sep=' ' quality_summaries} > quality_summary.tsv
      cat ${sep=' ' completeness} > completeness.tsv

      tar -xvf ${workflowscripts}
      Rscript --vanilla --slave Rscripts/parse_checkv_results.R quality_summary.tsv completeness.tsv vambviruses
  }
  output{
    File checkv_HQ = "vambviruses.HQ.tsv"
    File checkv_MQ = "vambviruses.MQ.tsv"
    File checkv_LQ = "vambviruses.LQ.tsv"

  }
  runtime{
      docker: "gcr.io/microbiome-xavier/r-docker:010319"
      disks: "local-disk 250 HDD"
      cpu : "1"
      preemptible: 2
  }
}


task pack_viruses{
   input{
    Array[File] viralbins
    File checkv_HQ
    File checkv_MQ
    File checkv_LQ
  }
  command <<<
      cat ~{sep=' ' viralbins} > all_contigs.fna.gz
      pigz -d all_contigs.fna.gz

      ### Select HQ and MQ viruses
      cut -f1 ~{checkv_HQ} > HQbins
      cut -f1 ~{checkv_MQ} > MQbins
      cut -f1 ~{checkv_LQ} > LQbins 

      ### 
      seqtk subseq all_contigs.fna HQbins > HQbins.fna
      seqtk subseq all_contigs.fna MQbins > MQbins.fna
      seqtk subseq all_contigs.fna LQbins > LQbins.fna
      
      ### Combine results 
      cat HQdir_dereplicated/* > HQbins.nr.fna 
      cat MQdir_dereplicated/* > MQbins.nr.fna
      pigz HQbins.nr.fna MQbins.nr.fna LQbins.fna MQbins.fna HQbins.fna
      cat ~{checkv_HQ} ~{checkv_MQ} ~{checkv_LQ} > checkv_quality_summary.tsv

      mkdir checkv
      mv *fna.gz checkv
      mv checkv_quality_summary.tsv checkv
      tar -zcvf checkvresults.tar.gz checkv
  >>>
  output{
      File Checkvresults = "checkvresults.tar.gz"
  }

  runtime{
      docker: "biojokke/coverm:0.6.1"
      disks: "local-disk 250 HDD"
      memory: "40GB"
      cpu : "32"
  }
}