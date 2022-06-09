version development

workflow viral_evaluation{

  input{
  ### Inputs to workflow 
  File viralreps
  File mgvrepo
  File mgv_proteins
  File workflowscripts
  File VOG77
  }

  call virfinder{
          input:
            viralreps = viralreps
  }

  ### 
  call prodigal_and_scan{
          input:
          viralreps = viralreps,
          mgvrepo = mgvrepo,
          Virfinder = virfinder.Virfinder
  }

  call viral_marker_tree{
          input:
          viralevaluation = prodigal_and_scan.viralevaluation,
          mgv_proteins = mgv_proteins,
          workflowscripts = workflowscripts,
          VOG77 = VOG77
  }

  call aai_clustering{
          input:
          viralevaluation = prodigal_and_scan.viralevaluation,
          mgv_proteins = mgv_proteins,
          workflowscripts = workflowscripts
  }
  
  output {
    File viruseval = prodigal_and_scan.viralevaluation 
  }
  
}


task virfinder{
    input{
      File viralreps
    }
    command <<<     
    tar -xvf ~{viralreps}

    Rscript /usr/local/bin/run_virfinder.Rscript ProVirusresults/viralreps.allviruses.fna viralreps.virfinder.tsv
    >>>
    output{
    File Virfinder = "viralreps.virfinder.tsv"
    }
    runtime{
    docker: "quay.io/fhcrc-microbiome/virfinder:v1.1--0"
    disks: "local-disk 250 HDD"
    memory: "40GB"
    cpu : "1"
    preemptible: 1
  }

}

task prodigal_and_scan{
  input{
    File viralreps
    File mgvrepo
    File Virfinder
  }
  command{
    tar -xvf ${mgvrepo}
    cd viral_detection_pipeline
    mv ${Virfinder} output/virfinder.tsv
    tar -xvf ${viralreps}
    mv ProVirusresults/viralreps.allviruses.fna input/viralreps.fna

    wget -O input/imgvr.hmm.gz https://img.jgi.doe.gov//docs/final_list.hmms.gz
    wget -O input/pfam.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz
    gunzip input/imgvr.hmm.gz
    gunzip input/pfam.hmm.gz

    prodigal -i input/viralreps.fna -a input/viralreps.faa -d input/viralreps.ffn -p meta -f gff > input/viralreps.gff
    hmmsearch -Z 1 --cpu 32 --noali --tblout output/imgvr.out input/imgvr.hmm input/viralreps.faa 
    hmmsearch -Z 1 --cut_tc --cpu 32 --noali --tblout output/pfam.out input/pfam.hmm input/viralreps.faa

    python count_hmm_hits.py input/viralreps.fna input/viralreps.faa output/imgvr.out output/pfam.out > output/hmm_hits.tsv
    python strand_switch.py input/viralreps.fna input/viralreps.faa > output/strand_switch.tsv
    python master_table.py output/hmm_hits.tsv output/virfinder.tsv output/strand_switch.tsv > output/master_table.tsv
    python viral_classify.py --features output/master_table.tsv --in_base input/viralreps --out_base output/viral
    
    gzip input/viralreps.faa input/viralreps.ffn input/viralreps.gff input/viralreps.fna output/imgvr.out output/pfam.out 
    rm input/imgvr.hmm 
    rm input/pfam.hmm

    ### Pack results
    cd ..
    tar -czvf MGVeval.tar.gz viral_detection_pipeline 
  }
  output{
    File viralevaluation = "MGVeval.tar.gz" 
  }

  runtime{
      docker: "quay.io/joacjo/vamb_virus:latest"
      disks: "local-disk 250 HDD"
      memory: "60GB"
      cpu : "40"
      preemptible: 1
  }
}



task viral_marker_tree{
  input{
    File viralevaluation
    File mgv_proteins
    File workflowscripts
    File VOG77
  }
  command{
    tar -xvf ${workflowscripts}
    tar -xvf ${viralevaluation} 
    gunzip -c ${mgv_proteins} > mgv_proteins.faa
    gunzip -c viral_detection_pipeline/input/viralreps.faa.gz > viralreps.faa
    cat viralreps.faa mgv_proteins.faa > viral.mgv.faa

    cp ${VOG77} VOG77.hmm
    python pyscripts/ViralMarkerTree77.py --in_faa viral.mgv.faa --out_dir markerout --threads 39
    
    ### Pack results
    tar -czvf viralmarker.tar.gz markerout
  }
  output{
    File viralmarker = "viralmarker.tar.gz" 
  }

  runtime{
      docker: "quay.io/joacjo/marker_tree:latest"
      disks: "local-disk 250 HDD"
      memory: "80GB"
      cpu : "40"
      preemptible: 1
  }
}


task aai_clustering{
  input{
    File viralevaluation
    File mgv_proteins
    File workflowscripts
  }
  command{
    tar -xvf ${workflowscripts}
    tar -xvf ${viralevaluation} 
    gunzip -c ${mgv_proteins} > mgv_proteins.faa
    gunzip -c viral_detection_pipeline/input/viralreps.faa.gz > viralreps.faa
    cat viralreps.faa mgv_proteins.faa > viral.mgv.faa

    diamond makedb --in viral.mgv.faa --db viral_proteins --threads 40
    diamond blastp --query viral.mgv.faa --db viral_proteins --out blastp.tsv \
    --outfmt 6 --evalue 1e-5 --max-target-seqs 10000 --query-cover 50 --subject-cover 50 --threads 40
    python pyscripts/calculate_genome_AAI.py --in_faa viral.mgv.faa --in_blast blastp.tsv --out_tsv aai.tsv

    python pyscripts/filter_aai.py --in_aai aai.tsv --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv genus_edges.tsv
    python pyscripts/filter_aai.py --in_aai aai.tsv --min_percent_shared 10 --min_num_shared 8 --min_aai 20 --out_tsv family_edges.tsv

    mcl genus_edges.tsv -te 8 -I 2.0 --abc -o genus_clusters.tsv
    mcl family_edges.tsv -te 8 -I 1.2 --abc -o family_clusters.tsv

    gzip blastp.tsv 
    gzip aai.tsv

    mkdir aai_clustering
    mv *tsv aai_clustering
    mv *tsv.gz aai_clustering

    ### Pack results
    tar -czvf aai.clusteringtar.gz aai_clustering

  }
  output{
    File genusclustering = "aai.clusteringtar.gz" 
  }

  runtime{
      docker: "quay.io/joacjo/marker_tree:latest"
      disks: "local-disk 250 HDD"
      memory: "120GB"
      cpu : "40"
      preemptible: 1
  }
}

