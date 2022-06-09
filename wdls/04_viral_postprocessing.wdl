  version development

  workflow virusMAGpostprocessing{

    input{
    ### Inputs to workflow 
    File checkvresults
    File checkmresults
    File msps
    File vamboutput
    File TPMmatrix
    File Genecatalogue
    File Proteincatalogue
    File workflowscripts
    File MGVdb
    File MGVdbFai
    File SpacerDB
    File PHROGDB
    File PFAMDB
    File MAGtax
    File MAGcrispr
    }

    call derepblast_viruses{
            input:
            checkvresults = checkvresults,
            Genecatalogue = Genecatalogue,
            workflowscripts = workflowscripts,
            MGVdb = MGVdb
    }

    call crisprhost{
            input:
            SpacerDB = SpacerDB,
            virusproducts = derepblast_viruses.virusproducts,
            workflowscripts = workflowscripts
    }

    call crisprhost_MAG{
            input:
            spacerhits = crisprhost.crisprblast,
            SpacerDB = SpacerDB,
            MAGcrispr = MAGcrispr,
            virusproducts = derepblast_viruses.virusproducts,
            workflowscripts = workflowscripts,
            MAGtax = MAGtax,
    }

    call virus_MGV_clustering{
            input:
              virusproducts = derepblast_viruses.virusproducts,
              workflowscripts = workflowscripts,
              MGVdb = MGVdb
    }

    call MGV_genecatalogue_coverage{
            input:
              Genecatalogue = Genecatalogue,
              workflowscripts = workflowscripts,
              MGVdb = MGVdb,
              MGVdbFai = MGVdbFai
    }

    call virfinder{
            input:
              virusproducts = derepblast_viruses.virusproducts,
              workflowscripts = workflowscripts,
    }

    call getvirusTPM{
            input:
            virusproducts = derepblast_viruses.virusproducts,
            virfinder = virfinder.Virfinder,
            TPMmatrix = TPMmatrix,
            workflowscripts = workflowscripts
      }

     call MAG_virus_prophages{
             input:
             virusreps = derepblast_viruses.virusreps,
             vamboutput = vamboutput,
             allcheckm = checkmresults,
             workflowscripts = workflowscripts
       }

     call viral_protein_annotation{
            input:
            Proteincatalogue = Proteincatalogue,
            workflowscripts = workflowscripts,
            PHROGDB = PHROGDB
     }

     call pfam_COG_annotation{
            input:
            Proteincatalogue = Proteincatalogue,
            workflowscripts = workflowscripts,
            PFAMDB = PFAMDB
     }
     call MAGMSPs{
             input:
             msps = msps,
             vamboutput = vamboutput,
             allcheckm = checkmresults,
             Genecatalogue = Genecatalogue,  
             workflowscripts = workflowscripts,
             MAGtax = MAGtax
     }
    
    output {
      File VirusResults = getvirusTPM.VirusResults
      File VirusProphages = MAG_virus_prophages.prophages
      File VirusHostMAG = crisprhost_MAG.hostassignment
      File MSPMAGs = MAGMSPs.MSPMAGs
    }
  }


task derepblast_viruses{
  input{
    File checkvresults
    File Genecatalogue
    File workflowscripts
    File MGVdb
  }
  command{
      ### 
      tar -xvf ${checkvresults}
      cat checkv/HQbins.fna.gz checkv/MQbins.fna.gz checkv/LQbins.fna.gz > HQMQLQ.bins.fna.gz
      pigz -d HQMQLQ.bins.fna.gz
      ls /app

      ### 
      makeblastdb -in HQMQLQ.bins.fna -dbtype nucl -out HQMQLQdb
      blastn -query HQMQLQ.bins.fna -db HQMQLQdb -outfmt '6 std qlen slen' \
      -max_target_seqs 10000 -perc_identity 90 \
      -out virusxvirus.m6 -num_threads 40

      python /app/anicalc.py -i virusxvirus.m6 -o virusxvirus.ani.tsv
      python /app/aniclust.py --fna HQMQLQ.bins.fna --ani virusxvirus.ani.tsv --out virus.clusters --min_ani 95 --min_qcov 0 --min_tcov 85
      
      cut -f1 virus.clusters  > viralreps.txt
      seqtk subseq HQMQLQ.bins.fna viralreps.txt > viralreps.fna 

      ### Combine with MGVdb 
      pigz -fdc ${MGVdb} > mgvreps.fna

      #cat viralreps.fna mgvreps.fna > viruses.fna
      cat viralreps.fna > viruses.fna
      makeblastdb -in viruses.fna -dbtype nucl -out HQMQLQdbreps 

      ### Blast Genecatalogue against the viral references
      gunzip -fdc ${Genecatalogue} > MGX.genecat.fna

      blastn -task megablast -db HQMQLQdbreps -query MGX.genecat.fna -out HQMQLQ.genecat.m6 -outfmt 6 -num_threads 40 \
      -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1

      tar -czvf virus.clustering.tar.gz viralreps.fna virus.clusters HQMQLQ.genecat.m6
  }
  output{
    File virusreps = "viralreps.fna"
    File virusclusters = "virus.clusters"
    File blastresults = "HQMQLQ.genecat.m6"
    File virusproducts = "virus.clustering.tar.gz"
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
        File virfinder
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
        Rscript --vanilla --slave Rscripts/get_viral_abundance_matrix.R HQMQLQ.genecat.m6 TPMmatrix.txt matrices

        ### pack results
        mkdir virusresults
        mv matrices/* virusresults
        mv viralreps.fna virusresults 
        mv virus.clusters virusresults
        mv HQMQLQ.genecat.m6 virusresults
        mv ${virfinder} virusresults
        tar -czvf virus.matrices.tar.gz virusresults
    }
    output{
    File TPM = "virusresults/combined_viralTPM.txt"
    File Frequency = "virusresults/combined_viralFrequency.txt"
    File VirusResults = "virus.matrices.tar.gz"
    }
    runtime{
      docker: "rocker/tidyverse:latest"
      disks: "local-disk 250 HDD"
      memory: "60GB"
      cpu : "4"
      preemptible: 1
    }
}


task MAGMSPs{
    input{
      File msps
      File vamboutput
      File allcheckm
      File Genecatalogue
      File workflowscripts
      File MAGtax
    }
    command <<<
        tar -xvf ~{vamboutput}
        tar -xvf ~{workflowscripts}
        ls vamb/bins/*fna.gz

        ### Assuming that bins are splitted already by sample and located under vamb/bins
        grep "::" ~{allcheckm} > allcheckm.txt
        awk '$13>=50 && $14<=10' allcheckm.txt > NCMQMAGS.txt
        
        mkdir tmp tmp2
        for MAG in `awk '{print $1}' NCMQMAGS.txt`; do 
            cp vamb/bins/${MAG}.fna.gz tmp
        done

        ### Concatenate, index and makeblastDB of MAGs
        ls tmp | head 
        cat tmp/*.fna.gz > NCMQMAGS.fna.gz 
        pigz -d NCMQMAGS.fna.gz
        makeblastdb -in NCMQMAGS.fna -dbtype nucl -out tmp2/MAGdb

        pigz -fdc ~{Genecatalogue} > genecat.fna
        mv ~{msps} combined_msp_gene_table.txt 
        mkdir msps
        ### Blast core-genes of each MSP against the MAGdb 
        for msp in `cut -f1 combined_msp_gene_table.txt | sort -u`; do 
            grep ${msp} combined_msp_gene_table.txt | grep -w core | cut -f4 > msps/${msp}.core_genes.txt
            seqtk subseq genecat.fna msps/${msp}.core_genes.txt > msps/${msp}.core_genes.fna
            blastn -task megablast -db tmp2/MAGdb -query msps/${msp}.core_genes.fna -out msps/${msp}.core_genes.MAG.m6 -outfmt '6 std qlen slen' \
            -num_threads 24 -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1
        done
        python pyscripts/MSP_MAG.py -i tmp -m msps -o MSP_MAGs.txt -g ~{MAGtax}

        tar -czvf MSP_MAG.tar.gz MSP_MAGs.txt msps

    >>>
    output{
    File MSPMAGs = "MSP_MAG.tar.gz"
    }
    runtime{
    docker: "quay.io/joacjo/vamb_virus"
    disks: "local-disk 250 HDD"
    memory: "40GB"
    cpu : "24"
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
      pigz -fdc ~{MGVdb} > mgvreps.fna
      cat viralreps.fna mgvreps.fna > viruses.fna

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
      tar -czvf MGVclustering.tar.gz MGVclustering
    >>>
    output{
    File MGVcluster = "MGVclustering.tar.gz"
    }
    runtime{
    docker: "quay.io/joacjo/vamb_virus"
    disks: "local-disk 250 HDD"
    memory: "80GB"
    cpu : "40"
    preemptible: 1
  }

}

task MGV_genecatalogue_coverage{
    input{
      File Genecatalogue
      File workflowscripts
      File MGVdb
      File MGVdbFai
    }
    command{

        pigz -fdc ${MGVdb} > mgvreps.fna
        makeblastdb -in mgvreps.fna -dbtype nucl -out MGVreps 

        ### Blast Genecatalogue against the viral references
        blastn -task megablast -db MGVreps -query ${Genecatalogue} -out MGV.genecat.m6 -outfmt 6 -num_threads 40 \
        -evalue 0.01 -qcov_hsp_perc 80 -perc_identity 90 -word_size 20 -reward 2 -penalty -3 -max_target_seqs 10 -max_hsps 1
        sort -k2 MGV.genecat.m6 > MGV.genecat.sorted.m6

        python /app/genecat_genome_coverage.py -i MGV.genecat.sorted.m6 -f ${MGVdbFai} -o MGV.genome.cov.txt -l 100
        tar -czvf MGV.coverage.tar.gz MGV.genecat.sorted.m6 MGV.genome.cov.txt
    }
    output{
      File MGVgenecat = "MGV.coverage.tar.gz"
    }
    runtime{
        docker: "quay.io/joacjo/vamb_virus:latest"
        disks: "local-disk 250 HDD"
        memory: "60GB"
        cpu : "40"
        preemptible: 1
    }
}


task virfinder{
    input{
      File virusproducts
      File workflowscripts
    }
    command <<<
      tar -xvf ~{virusproducts}
      tar -xvf ~{workflowscripts}
      
      ### Combine with MGVdb 
      head viralreps.fna
      Rscript /usr/local/bin/run_virfinder.Rscript viralreps.fna viralreps.virfinder.tsv

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

task viral_protein_annotation{
    input{
      File Proteincatalogue
      File workflowscripts
      File PHROGDB
    }
    command <<<
      tar -xvf ~{workflowscripts}
      tar -xvf ~{PHROGDB}

      ### Run MMSeqs2 search against the PHROG database
      cd phrogs_mmseqs_db
      gunzip -fdc ~{Proteincatalogue} > MGX.genecat.faa
      mmseqs createdb MGX.genecat.faa target_seq
      
      mmseqs search phrogs_profile_db target_seq results_mmseqs ./tmp -s 5.7
      mmseqs createtsv phrogs_profile_db target_seq results_mmseqs MGX.genecat.phrogs.tsv
      
      ls 
      #tar -czvf viral.protein.annotation.tar.gz MGX.genecat.phrogs.tsv
    >>>
    output{
    File phrogs = "phrogs_mmseqs_db/MGX.genecat.phrogs.tsv"
    }
    runtime{
    docker: "soedinglab/mmseqs2"
    disks: "local-disk 250 HDD"
    memory: "80GB"
    cpu : "40"
    preemptible: 1
  }
}


task crisprhost{
      input{
        File SpacerDB
        File virusproducts
        File workflowscripts
      }
      command <<<
      tar -xvf ~{SpacerDB}
      tar -xvf ~{virusproducts}
      tar -xvf ~{workflowscripts}
      mkdir -p tmp 
      makeblastdb -in viralreps.fna -dbtype nucl -out tmp/viraldb
      grep '>' viralreps.fna > viralreps.txt

      gunzip SpacersDB.fasta.gz 
      blastn -task blastn-short -num_threads 40 \
      -evalue 0.01 \
      -db tmp/viraldb \
      -query SpacersDB.fasta \
      -out spacerhits.m6 \
      -outfmt '6 std qlen slen'

      ls

      >>>
      output{
      File crisprblast = "spacerhits.m6"
      }
      runtime{
      docker: "quay.io/joacjo/vamb_virus"
      disks: "local-disk 250 HDD"
      memory: "40GB"
      cpu : "40"
      preemptible: 1
    }
  }

task crisprhost_MAG{
      input{
        File spacerhits
        File SpacerDB
        File MAGcrispr
        File virusproducts
        File workflowscripts
        File MAGtax
      }
      command <<<
      tar -xvf ~{SpacerDB}
      tar -xvf ~{MAGcrispr} # all.spacers.fna all.crtspacers.fna all.crisprs.tab
      tar -xvf ~{virusproducts}
      tar -xvf ~{workflowscripts}
      mkdir -p tmp 
      makeblastdb -in viralreps.fna -dbtype nucl -out tmp/viraldb
      grep '>' viralreps.fna > viralreps.txt

      cat all.spacers.fna all.crtspacers.fna  > mag.spacers.fna
      blastn -task blastn-short -num_threads 40 \
      -evalue 0.01 \
      -db tmp/viraldb \
      -query mag.spacers.fna \
      -out MAGspacerhits.m6 \
      -outfmt '6 std qlen slen'

      ls
      python pyscripts/crispr_host_annotation.py -g CrisprOpenDB.tax.tsv -m ~{MAGtax} -v viralreps.txt -o virus_MAG_crispr.txt -b MAGspacerhits.m6 -c ~{spacerhits}

      >>>
      output{
      File hostassignment = "virus_MAG_crispr.txt"
      File crisprblast = "MAGspacerhits.m6"
      }
      runtime{
      docker: "quay.io/joacjo/vamb_virus"
      disks: "local-disk 250 HDD"
      memory: "40GB"
      cpu : "40"
      preemptible: 1
    }
  }

task pfam_COG_annotation{
    input{
      File Proteincatalogue
      File workflowscripts
      File PFAMDB
    }
    command <<<
      tar -xvf ~{workflowscripts}
      tar -xvf ~{PFAMDB}

      ### Run MMSeqs2 search against the PHROG database
      gunzip -fdc ~{Proteincatalogue} > MGX.genecat.faa
      hmmsearch -E 0.01 --tblout MGX.genecat.pfam.tbl --cpu 40 Pfam-A.hmm MGX.genecat.faa > stdout.txt
    >>>
    output{
    File pfamtbl = "MGX.genecat.pfam.tbl"
    }
    runtime{
    docker: "staphb/hmmer"
    disks: "local-disk 250 HDD"
    memory: "80GB"
    cpu : "40"
    preemptible: 1
  }
}